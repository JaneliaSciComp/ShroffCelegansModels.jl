"""
Scans every dataset listed in the active config_path and records mtimes of two
kinds of inputs to the recompute pipeline:

- "annotation" — `Decon_reg_*/.../integrated_annotation/annotations.csv`
- "lattice"    — `Decon_reg_*/.../lattice_final/lattice.csv` plus the
                 `model_crossSections/latticeCrossSection_*.csv` files

Writes two HDF5 files into the directory given by MODIFIED_TIMES_DIR
(default: /data/annotations/modified_times):

  - modified_times_YYYY_mm_dd_HHMMSS.h5  (timestamped snapshot, kept as history)
  - modified_times.h5                    (overwritten each run, "latest" pointer)

Both files use a `kind/group/dataset_index` HDF5 layout. The prior layout was
flat `group/dataset_index` for annotation mtimes only; that legacy shape is
detected on load and treated as "no prior data" so the first new snapshot
self-heals without manual PVC cleanup.

After writing, compares the new mtimes against the prior `modified_times.h5`.
If any (group, dataset, timepoint) mtime advanced — for either kind — drops a
JSON marker at `MODIFIED_TIMES_DIR/pending_recompute` carrying the change
count, a sample, and the sorted unique kinds that changed. The marker pattern
avoids the need for in-cluster K8s API dispatch (which would require RBAC the
namespace user can't grant).

A separate frequently-running CronJob (`run_recompute_if_needed.jl`) picks up
the marker and runs the recompute pipeline.
"""

using Dates: now, format
using HDF5: HDF5, h5open
using ShroffCelegansModels
using ShroffCelegansModels.JSON3
using ShroffCelegansModels.MIPAVIO: get_annotation_modified_times_unix,
                                    get_lattice_modified_times_unix,
                                    save_all_modified_times_unix

using ShroffCelegansModels: read_config_json

# A single (group, dataset_index, timepoint_index) tuple whose mtime advanced
# between two snapshots, tagged by `kind` (:annotation or :lattice).
# timepoint_index is a 1-based offset within the dataset's stored mtimes vector.
struct MtimeChange
    kind::Symbol
    group::String
    dataset_index::Int
    timepoint_index::Int
    old::Float64
    new::Float64
end

# Read the per-kind mtimes from a snapshot file into the shape
# `Dict{group => Vector{Vector{Float64}}}` (parallel to what
# `get_*_modified_times_unix(datasets)` returns).
#
# Returns `nothing` if the file is absent, has no top-level group for `kind`,
# OR if it's in the legacy flat layout (group keys point directly to dataset
# groups instead of a kind group). Legacy detection: any top-level child whose
# own children have integer-parseable names is a dataset group, not a kind
# group → legacy.
function read_prior_mtimes(path::AbstractString, kind::AbstractString)::Union{Nothing, Dict{String, Vector{Vector{Float64}}}}
    isfile(path) || return nothing
    h5open(path, "r") do f
        top = collect(keys(f))
        isempty(top) && return nothing
        # Legacy detection: take any top-level group and check whether its
        # children parse as integers (= dataset indices). If so, the file
        # predates the kind layout and we have no prior data for either kind.
        first_top = f[first(top)]
        if first_top isa HDF5.Group && !isempty(keys(first_top))
            sample = first(keys(first_top))
            if tryparse(Int, sample) !== nothing
                return nothing
            end
        end
        haskey(f, kind) || return nothing
        kind_group = f[kind]
        out = Dict{String, Vector{Vector{Float64}}}()
        for group_name in keys(kind_group)
            g = kind_group[group_name]
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
    kind::Symbol,
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
                    push!(changes, MtimeChange(kind, group, ds_idx, tp, o, n))
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
    annotation_mtimes = get_annotation_modified_times_unix(datasets)
    @info "Scanning lattice mtimes"
    lattice_mtimes = get_lattice_modified_times_unix(datasets)

    timestamp = format(now(), "yyyy_mm_dd_HHMMSS")
    snapshot_path = joinpath(output_dir, "modified_times_$(timestamp).h5")
    latest_path = joinpath(output_dir, "modified_times.h5")

    @info "Reading prior modified_times for diff" latest_path
    prior_annotation = read_prior_mtimes(latest_path, "annotation")
    prior_lattice    = read_prior_mtimes(latest_path, "lattice")
    annotation_changes = prior_annotation === nothing ? MtimeChange[] :
                         diff_mtimes(:annotation, annotation_mtimes, prior_annotation)
    lattice_changes    = prior_lattice === nothing ? MtimeChange[] :
                         diff_mtimes(:lattice, lattice_mtimes, prior_lattice)
    changes = vcat(annotation_changes, lattice_changes)

    if prior_annotation === nothing && prior_lattice === nothing
        @info "No prior snapshot for either kind — first run / post-migration, no diff to compute"
    else
        @info "Diff complete" annotation_changes=length(annotation_changes) lattice_changes=length(lattice_changes) examples=first(changes, min(5, length(changes)))
    end

    @info "Writing snapshot" snapshot_path
    write_atomic(snapshot_path) do tmp
        save_all_modified_times_unix(datasets; filepath=tmp, annotation_mtimes, lattice_mtimes)
    end

    @info "Updating latest" latest_path
    write_atomic(latest_path) do tmp
        save_all_modified_times_unix(datasets; filepath=tmp, annotation_mtimes, lattice_mtimes)
    end

    if !isempty(changes)
        marker_path = joinpath(output_dir, "pending_recompute")
        kinds = sort(unique(string(c.kind) for c in changes))
        marker = Dict(
            "triggered_at" => format(now(), "yyyy-mm-ddTHH:MM:SS"),
            "change_count" => length(changes),
            "kinds" => kinds,
            "examples" => [
                Dict(
                    "kind" => string(c.kind),
                    "group" => c.group,
                    "dataset_index" => c.dataset_index,
                    "timepoint_index" => c.timepoint_index,
                    "old_mtime" => isnan(c.old) ? nothing : c.old,
                    "new_mtime" => c.new,
                )
                for c in first(changes, min(20, length(changes)))
            ],
        )
        @info "Writing recompute marker" marker_path change_count=length(changes) kinds
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
