"""
Scans every dataset listed in the active config_path and records the mtime of
each Decon_reg_*/.../integrated_annotation/annotations.csv file. Used by the
OpenShift CronJob to track when integrated annotations were last edited.

Writes two HDF5 files into the directory given by MODIFIED_TIMES_DIR
(default: /data/annotations/modified_times):

  - modified_times_YYYY_mm_dd_HHMMSS.h5  (timestamped snapshot, kept as history)
  - modified_times.h5                    (overwritten each run, "latest" pointer)

After writing, compares the new mtimes against the prior `modified_times.h5`.
If any (group, dataset, timepoint) mtime advanced, drops a JSON marker at
`MODIFIED_TIMES_DIR/pending_recompute`. A separate frequently-running CronJob
(`run_recompute_if_needed.jl`) picks up the marker and runs the recompute
pipeline. The marker pattern avoids the need for in-cluster K8s API dispatch
(which would require RBAC the namespace user can't grant).
"""

using Dates: now, format
using HDF5: h5open
using ShroffCelegansModels
using ShroffCelegansModels.JSON3
using ShroffCelegansModels.MIPAVIO: get_modified_times_unix, save_modified_times_unix

include(joinpath(@__DIR__, "..", "..", "src", "demo_averaging", "read_config_json.jl"))

# A single (group, dataset_index, timepoint_index) tuple whose mtime advanced
# between two snapshots. timepoint_index is 1-based offset within the dataset's
# stored mtimes vector (matches the existing modified_times.h5 layout).
struct MtimeChange
    group::String
    dataset_index::Int
    timepoint_index::Int
    old::Float64
    new::Float64
end

# Read a prior modified_times.h5 into Dict{group => Vector{Vector{Float64}}},
# parallel to the in-memory shape `get_modified_times_unix` returns. Returns
# `nothing` if the file doesn't exist (first-ever run).
function read_prior_mtimes(path::AbstractString)::Union{Nothing, Dict{String, Vector{Vector{Float64}}}}
    isfile(path) || return nothing
    h5open(path, "r") do f
        out = Dict{String, Vector{Vector{Float64}}}()
        for group_name in keys(f)
            g = f[group_name]
            indices = sort(parse.(Int, collect(keys(g))))
            vecs = Vector{Vector{Float64}}(undef, length(indices))
            for (slot, idx) in enumerate(indices)
                vecs[slot] = Float64.(read(g[string(idx)]))
            end
            out[group_name] = vecs
        end
        return out
    end
end

# Per-element comparison: a change requires the new mtime to be finite AND
# strictly greater than the old (treating prior NaN/absent as "no record").
# Datasets/timepoints present only in `new` are reported as changes from NaN.
function diff_mtimes(
    new::Dict{String, Vector{Vector{Float64}}},
    old::Union{Nothing, Dict{String, Vector{Vector{Float64}}}},
)::Vector{MtimeChange}
    changes = MtimeChange[]
    for (group, new_vecs) in new
        old_vecs = old === nothing ? nothing : get(old, group, nothing)
        for (ds_idx, new_v) in enumerate(new_vecs)
            old_v = old_vecs === nothing || ds_idx > length(old_vecs) ? nothing : old_vecs[ds_idx]
            for tp in eachindex(new_v)
                n = new_v[tp]
                isnan(n) && continue
                o = old_v === nothing || tp > length(old_v) ? NaN : old_v[tp]
                if isnan(o) || n > o
                    push!(changes, MtimeChange(group, ds_idx, tp, o, n))
                end
            end
        end
    end
    return changes
end

# Write `dst` atomically: invoke `f(tmp_path)` to produce the file, then `mv` it
# into place. POSIX rename() is atomic on the same filesystem, so readers always
# see either the previous version or the new one — never a partial write.
function write_atomic(f, dst::AbstractString)
    tmp = string(dst, ".tmp.", getpid(), ".", time_ns())
    try
        f(tmp)
        mv(tmp, dst; force=true)
    catch
        isfile(tmp) && rm(tmp; force=true)
        rethrow()
    end
end

function main()
    output_dir = get(ENV, "MODIFIED_TIMES_DIR", "/data/annotations/modified_times")
    mkpath(output_dir)

    config_path = ShroffCelegansModels.config_path
    @info "Loading datasets from config" config_path
    _, _, datasets = read_config_json(config_path)
    @info "Loaded datasets" groups=length(datasets) total=sum(length, values(datasets))

    @info "Scanning integrated_annotation mtimes"
    modified_times = get_modified_times_unix(datasets)

    timestamp = format(now(), "yyyy_mm_dd_HHMMSS")
    snapshot_path = joinpath(output_dir, "modified_times_$(timestamp).h5")
    latest_path = joinpath(output_dir, "modified_times.h5")

    @info "Reading prior modified_times for diff" latest_path
    prior = read_prior_mtimes(latest_path)
    changes = prior === nothing ? MtimeChange[] : diff_mtimes(modified_times, prior)
    if prior === nothing
        @info "No prior snapshot found — first run, no diff to compute"
    else
        @info "Diff complete" changes=length(changes) examples=first(changes, min(5, length(changes)))
    end

    @info "Writing snapshot" snapshot_path
    write_atomic(snapshot_path) do tmp
        save_modified_times_unix(datasets, modified_times; filepath=tmp)
    end

    @info "Updating latest" latest_path
    write_atomic(latest_path) do tmp
        save_modified_times_unix(datasets, modified_times; filepath=tmp)
    end

    if !isempty(changes)
        marker_path = joinpath(output_dir, "pending_recompute")
        marker = Dict(
            "triggered_at" => format(now(), "yyyy-mm-ddTHH:MM:SS"),
            "change_count" => length(changes),
            "examples" => [
                Dict(
                    "group" => c.group,
                    "dataset_index" => c.dataset_index,
                    "timepoint_index" => c.timepoint_index,
                    "old_mtime" => isnan(c.old) ? nothing : c.old,
                    "new_mtime" => c.new,
                )
                for c in first(changes, min(20, length(changes)))
            ],
        )
        @info "Writing recompute marker" marker_path change_count=length(changes)
        write_atomic(marker_path) do tmp
            open(tmp, "w") do io
                JSON3.write(io, marker)
            end
        end
    end

    @info "Done"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
