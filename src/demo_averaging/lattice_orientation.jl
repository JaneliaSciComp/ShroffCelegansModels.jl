# Left-right (LR) orientation QC: detect `lattice.csv` files whose Left/Right
# seam-cell columns are swapped, by checking whether Cpaaaa (hyp7) sits on the
# expected side of the seam-cell plane. Reuses the untwisted annotation frame
# from `annotation_untwist.jl` (`untwist_annotation`) — no new axis geometry.
#
# Business logic (majority-vote reference sign, mismatch flags, HDF5 writing)
# lives in `web/scripts/check_lattice_orientation.jl`, not here — this file
# only extracts the raw per-timepoint signal, mirroring how `MIPAVIO.jl` keeps
# `get_lattice_modified_times_unix` separate from `save_modified_times.jl`'s
# diffing logic.

const CPAAAA_NAME_CANDIDATES = ("Cpaaaa", "hyp7_Cpaaaa", "hyp7")

# Resolve which key in a twisted-annotations dict corresponds to Cpaaaa/hyp7.
# The annotations-CSV `name` column may hold the lineage name or the MIPAV
# positional name, so try both directly before falling back to a
# case-insensitive scan of `cell_key.mapping`'s display names.
function resolve_cpaaaa_key(dict::AbstractDict, cell_key::CellKey)
    for c in CPAAAA_NAME_CANDIDATES
        haskey(dict, c) && return c
    end
    for (raw, display) in cell_key.mapping
        if any(c -> lowercase(display) == lowercase(c), CPAAAA_NAME_CANDIDATES)
            raw_str = string(raw)
            haskey(dict, raw_str) && return raw_str
        end
    end
    return nothing
end

# Sign of Cpaaaa's position along the DV axis of the untwisted annotation
# frame (see `untwist_annotation`), relative to the model's seam-cell plane.
# `NaN` for degenerate cases (exactly on the plane) or any failure.
function lattice_orientation_sign(model::AbstractCelegansModel, pt::Point)
    try
        u = untwist_annotation(model, pt)
        s = sign(u[2])
        return s == 0.0 ? NaN : Float64(s)
    catch
        return NaN
    end
end

# One dataset: a sign per timepoint (`NaN` for outliers / missing Cpaaaa
# annotation / build failures), mirroring
# `MIPAVIO.get_lattice_modified_times_unix(dataset)`.
function lattice_orientation_signs(dataset::NormalizedDataset)::Vector{Float64}
    mts = ModelTimeSeries(dataset)
    n = length(range(dataset.cell_key))
    map(1:n) do i
        try
            model = mts(i)
            ismissing(model) && return NaN
            dict = twisted_annotations(dataset, i)
            ismissing(dict) && return NaN
            key = resolve_cpaaaa_key(dict, dataset.cell_key)
            isnothing(key) && return NaN
            lattice_orientation_sign(model, dict[key])
        catch
            NaN
        end
    end
end

function lattice_orientation_signs(
    datasets::Dict{String, Vector{NormalizedDataset}}
)::Dict{String, Vector{Vector{Float64}}}
    Dict(k => lattice_orientation_signs.(v) for (k, v) in datasets)
end

# The resolved Cpaaaa/hyp7 annotation key for a dataset, taken from its first
# timepoint where it's resolvable (`nothing` if never found across the whole
# series). Purely a diagnostic value for the QC snapshot — doesn't need the
# lattice model, just the annotations CSV, so it's cheap to compute alongside
# `lattice_orientation_signs`.
function lattice_orientation_cpaaaa_key(dataset::NormalizedDataset)
    n = length(range(dataset.cell_key))
    for i in 1:n
        try
            dict = twisted_annotations(dataset, i)
            ismissing(dict) && continue
            key = resolve_cpaaaa_key(dict, dataset.cell_key)
            isnothing(key) || return key
        catch
        end
    end
    return nothing
end
