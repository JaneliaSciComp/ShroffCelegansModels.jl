"""
Scans every dataset listed in the active config_path and checks whether each
timepoint's tracked lattice is oriented correctly in the left-right (LR)
axis. The check exploits a fixed biological invariant: Cpaaaa (hyp7) sits
dorsal to the seam-cell plane in every correctly-tracked worm, so a
`lattice.csv` whose Left/Right seam-cell columns got swapped — for a whole
dataset, or for a single re-tracked timepoint — mirrors that sign.

There is no hardcoded "dorsal" ground truth: the reference sign is
self-calibrated as the majority sign across all datasets (one dataset, one
vote, so a single long time series can't dominate). Two independent flags
are computed against that reference:

- `matches_reference` (per dataset): the dataset's own majority sign agrees
  with the global reference — disagreement suggests the whole dataset has
  LR swapped.
- per-timepoint values that disagree with their *own* dataset's majority sign
  suggest a single mislabeled/re-tracked `lattice.csv`, independent of
  whether the dataset itself matches the global reference.

Writes `LATTICE_ORIENTATION_DIR/lattice_orientation.h5` (env var, default
/data/annotations/lattice_orientation) plus a timestamped snapshot copy, using
the same `kind/group/dataset_index` HDF5 shape as `modified_times.h5`
(`kind` = "lattice_orientation").
"""

using Dates: now, format
using HDF5: HDF5, h5open, create_group, attrs
using ShroffCelegansModels
using ShroffCelegansModels: read_config_json, NormalizedDataset,
                             lattice_orientation_signs, lattice_orientation_cpaaaa_key

# Majority value among {+1.0, -1.0}, ignoring NaN; NaN if nothing valid.
# Ties default to +1.0 — arbitrary but deterministic.
function majority_sign(xs)
    valid = filter(!isnan, xs)
    isempty(valid) && return NaN
    pos = count(==(1.0), valid)
    neg = count(==(-1.0), valid)
    return pos >= neg ? 1.0 : -1.0
end

# Write `dst` atomically: invoke `f(tmp_path)` to produce the file, then `mv`
# it into place, matching `save_modified_times.jl`'s convention.
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

function save_lattice_orientation(
    filepath::String,
    datasets::Dict{String, Vector{NormalizedDataset}},
    signs::Dict{String, Vector{Vector{Float64}}},
    representative::Dict{String, Vector{Float64}},
    cpaaaa_keys::Dict{String, Vector{Union{Nothing,String}}},
    reference_sign::Float64,
)
    h5open(filepath, "w") do h5f
        A = attrs(h5f)
        A["reference_sign"] = reference_sign
        A["generated_at"] = format(now(), "yyyy-mm-ddTHH:MM:SS")
        kind_group = create_group(h5f, "lattice_orientation")
        for group in keys(datasets)
            h5g = create_group(kind_group, group)
            for (k, ds) in enumerate(datasets[group])
                rep = representative[group][k]
                h5g[string(k)] = signs[group][k]
                dsattrs = attrs(h5g[string(k)])
                dsattrs["path"] = ds.path
                dsattrs["cell_key.name"] = ds.cell_key.name
                dsattrs["cell_key.start"] = ds.cell_key.start
                dsattrs["cell_key.end"] = ds.cell_key.stop
                dsattrs["cell_key.outliers"] = ds.cell_key.outliers
                dsattrs["cpaaaa_key"] = something(cpaaaa_keys[group][k], "")
                dsattrs["representative_sign"] = rep
                dsattrs["matches_reference"] = !isnan(rep) && rep == reference_sign
            end
        end
    end
end

function main()
    output_dir = get(ENV, "LATTICE_ORIENTATION_DIR", "/data/annotations/lattice_orientation")
    mkpath(output_dir)

    config_path = ShroffCelegansModels.config_path
    @info "Loading datasets from config" config_path
    _, _, datasets = read_config_json(config_path)
    @info "Loaded datasets" groups=length(datasets) total=sum(length, values(datasets))

    @info "Checking lattice LR orientation (Cpaaaa vs. seam plane)"
    signs = lattice_orientation_signs(datasets)
    cpaaaa_keys = Dict{String, Vector{Union{Nothing,String}}}(
        group => lattice_orientation_cpaaaa_key.(datasets[group]) for group in keys(datasets)
    )

    representative = Dict(group => majority_sign.(signs[group]) for group in keys(signs))
    all_representative = reduce(vcat, values(representative); init=Float64[])
    reference_sign = majority_sign(all_representative)
    mismatched = count(x -> !isnan(x) && x != reference_sign, all_representative)

    @info "Orientation check complete" reference_sign dataset_count=length(all_representative) mismatched_dataset_count=mismatched

    timestamp = format(now(), "yyyy_mm_dd_HHMMSS")
    snapshot_path = joinpath(output_dir, "lattice_orientation_$(timestamp).h5")
    latest_path = joinpath(output_dir, "lattice_orientation.h5")

    @info "Writing snapshot" snapshot_path
    write_atomic(snapshot_path) do tmp
        save_lattice_orientation(tmp, datasets, signs, representative, cpaaaa_keys, reference_sign)
    end

    @info "Updating latest" latest_path
    write_atomic(latest_path) do tmp
        save_lattice_orientation(tmp, datasets, signs, representative, cpaaaa_keys, reference_sign)
    end

    @info "Done"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
