# Left-right (LR) orientation QC: detect `lattice.csv` files whose Left/Right
# seam-cell columns are swapped, by checking whether Cpaaaa (hyp7) sits on the
# expected side of the seam-cell plane. Reuses the untwisted annotation frame
# from `annotation_untwist.jl` (`untwist_annotation`) — no new axis geometry.
#
# `untwist_annotation` returns Point3d(x, y, z) where `x` is the signed
# distance along the "right" transverse direction (LR axis) and `y` is the
# signed distance along the DV axis, both measured from the seam-cell plane
# at the nearest point along the central (AP) spline. Both are useful QC
# invariants: Cpaaaa (hyp7) is reliably dorsal (`y > 0`), and some other,
# strain-independent cells may instead be reliably lateralized (`x`
# consistently one sign) — either kind of signal can catch an LR swap, and
# the magnitude of either component is a confidence measure (a point sitting
# right on the plane, `magnitude ≈ 0`, is a weak/unreliable signal even if
# its sign happens to be "correct").
#
# Business logic (majority-vote reference sign, mismatch flags, HDF5 writing,
# candidate-cell ranking) lives in `web/scripts/check_lattice_orientation.jl`
# and `web/scripts/survey_cell_orientation_candidates.jl`, not here — this
# file only extracts the raw per-timepoint signal, mirroring how
# `MIPAVIO.jl` keeps `get_lattice_modified_times_unix` separate from
# `save_modified_times.jl`'s diffing logic.

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

# Sign+magnitude along both untwisted-frame axes for one annotation point,
# relative to the model's seam-cell plane at that AP position:
# - `lr_sign`/`lr_magnitude`: the LR (right-transverse, `x`) component.
# - `dv_sign`/`dv_magnitude`: the DV (normal, `y`) component.
# Signs are `NaN` for degenerate cases (exactly on an axis) or any failure;
# magnitudes are `NaN` only on failure (a point on the plane has magnitude 0
# with a NaN sign, which is exactly the "no signal" case callers should skip).
const NAN_AXES = (lr_sign = NaN, lr_magnitude = NaN, dv_sign = NaN, dv_magnitude = NaN)

function lattice_orientation_axes(model::AbstractCelegansModel, pt::Point)
    try
        u = untwist_annotation(model, pt)
        x, y = Float64(u[1]), Float64(u[2])
        lr_s, dv_s = sign(x), sign(y)
        return (
            lr_sign = lr_s == 0.0 ? NaN : lr_s,
            lr_magnitude = abs(x),
            dv_sign = dv_s == 0.0 ? NaN : dv_s,
            dv_magnitude = abs(y),
        )
    catch
        return NAN_AXES
    end
end

# Same as `lattice_orientation_axes`, batched over every cell in `dict` at
# one timepoint (a single `untwist_annotations` call instead of one per
# cell) — used by `dataset_cell_orientation_survey` to profile every
# annotated cell, not just Cpaaaa.
function lattice_orientation_axes_all(model::AbstractCelegansModel, dict::AbstractDict)
    names = collect(keys(dict))
    isempty(names) && return Dict{String, typeof(NAN_AXES)}()
    pts = [dict[n] for n in names]
    axes = try
        us = untwist_annotations(model, pts)
        map(us) do u
            x, y = Float64(u[1]), Float64(u[2])
            lr_s, dv_s = sign(x), sign(y)
            (
                lr_sign = lr_s == 0.0 ? NaN : lr_s,
                lr_magnitude = abs(x),
                dv_sign = dv_s == 0.0 ? NaN : dv_s,
                dv_magnitude = abs(y),
            )
        end
    catch
        fill(NAN_AXES, length(names))
    end
    return Dict(string(names[i]) => axes[i] for i in eachindex(names))
end

# One dataset, one cell (Cpaaaa/hyp7): per-timepoint axes (`NAN_AXES` for
# outliers / missing annotation / build failures), mirroring
# `MIPAVIO.get_lattice_modified_times_unix(dataset)`.
function lattice_orientation_series(dataset::NormalizedDataset)
    mts = ModelTimeSeries(dataset)
    n = length(range(dataset.cell_key))
    map(1:n) do i
        try
            model = mts(i)
            ismissing(model) && return NAN_AXES
            dict = twisted_annotations(dataset, i)
            ismissing(dict) && return NAN_AXES
            key = resolve_cpaaaa_key(dict, dataset.cell_key)
            isnothing(key) && return NAN_AXES
            lattice_orientation_axes(model, dict[key])
        catch
            NAN_AXES
        end
    end
end

function lattice_orientation_series(datasets::Dict{String, Vector{NormalizedDataset}})
    Dict(k => lattice_orientation_series.(v) for (k, v) in datasets)
end

# The resolved Cpaaaa/hyp7 annotation key for a dataset, taken from its first
# timepoint where it's resolvable (`nothing` if never found across the whole
# series). Purely a diagnostic value for the QC snapshot — doesn't need the
# lattice model, just the annotations CSV, so it's cheap to compute alongside
# `lattice_orientation_series`.
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

# Every cell's per-timepoint axes for one dataset — used to survey candidate
# anchor cells besides Cpaaaa (e.g. for strains where it isn't tracked).
# Timepoints where the model can't be built or the annotations file can't be
# read are skipped entirely (not NaN-filled), since they contribute no cell
# names at all; a cell simply has fewer observations for that dataset.
#
# The annotations-CSV `name` column is frequently a MIPAV positional/tracking
# placeholder (e.g. "C10", "C2") rather than a stable identity — the same
# placeholder resolves to entirely different cells in different datasets
# (confirmed on real data: "C10" maps to 4 unrelated cells, "C2" to 11).
# Unlike `resolve_cpaaaa_key` (which also accepts a few known direct-name
# candidates before falling back to `cell_key.mapping`, since Cpaaaa's exact
# identity is externally known), a raw name here is trusted only when this
# dataset's own `cell_key.mapping` explicitly resolves it — a name absent
# from the mapping is dropped rather than assumed to already be canonical,
# since that assumption can't be verified per-dataset and was previously
# found to admit placeholders whose mapping was simply incomplete.
function dataset_cell_orientation_survey(dataset::NormalizedDataset)
    mts = ModelTimeSeries(dataset)
    n = length(range(dataset.cell_key))
    mapping = dataset.cell_key.mapping
    result = Dict{String, Vector{typeof(NAN_AXES)}}()
    for i in 1:n
        try
            model = mts(i)
            ismissing(model) && continue
            dict = twisted_annotations(dataset, i)
            ismissing(dict) && continue
            per_cell = lattice_orientation_axes_all(model, dict)
            for (raw_name, ax) in per_cell
                sym = Symbol(raw_name)
                haskey(mapping, sym) || continue
                canonical = mapping[sym]
                push!(get!(() -> typeof(NAN_AXES)[], result, canonical), ax)
            end
        catch
        end
    end
    return result
end
