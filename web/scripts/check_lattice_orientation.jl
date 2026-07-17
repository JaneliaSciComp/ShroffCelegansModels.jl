"""
Scans every dataset listed in the active config_path and checks whether each
timepoint's tracked lattice is oriented correctly, using a set of named
"checks" — each a (cell name, axis) pair exploiting the fixed biological
invariant that a given cell sits reliably on one side (dorsal/ventral, or
left/right) of the seam-cell plane in every correctly-tracked worm. A
`lattice.csv` whose Left/Right seam-cell columns got swapped — for a whole
dataset, or for a single re-tracked timepoint — mirrors that sign.

`hyp7_Cpaaaa` (DV axis) is the original check. The other checks were found
by `web/scripts/survey_cell_orientation_candidates.jl`'s greedy set-cover: a
small set of cells, together resolvable across every strain, to serve as
fallback anchors where Cpaaaa isn't tracked (see CHECKS below).

There is no hardcoded "correct sign" ground truth: each check's reference
sign is self-calibrated as the majority sign across all datasets where that
check resolves (one dataset, one vote). Two independent flags are computed
against that reference, per check:

- `matches_reference` (per dataset): the dataset's own majority sign for
  this check agrees with the check's global reference — disagreement
  suggests the whole dataset has this axis swapped.
- per-timepoint values that disagree with their *own* dataset's majority
  sign for this check suggest a single mislabeled/re-tracked `lattice.csv`,
  independent of whether the dataset itself matches the global reference.
  The count of these is stored per (dataset, check) as
  `mismatched_timepoint_count`.

Cell resolution: `hyp7_Cpaaaa` uses `resolve_cpaaaa_key` (direct-name
candidates, falling back to `cell_key.mapping`); every other check uses
strict `cell_key.mapping`-only resolution (`resolve_mapped_key`), matching
how the survey validated them. A check is only written for a dataset when
it actually resolves there (representative_sign non-NaN) — most of the ~30
checks are irrelevant to any given dataset, so this keeps the snapshot and
the webapp focused on what's actually present. `link_eligible` records
whether the resolved name is one of `fix_annotation_ap_axis`'s selectable
cells for that dataset (`fix_ap_axis_selectable`), i.e. whether a deep link
into it with `?annotation=<name>` will work — always true for the
strictly-resolved checks, sometimes false for `hyp7_Cpaaaa` (found via
direct name in datasets whose `cell_key.mapping` doesn't cover it).

Writes `LATTICE_ORIENTATION_DIR/lattice_orientation.h5` (env var, default
/data/annotations/lattice_orientation) plus a timestamped snapshot copy.
"""

using Dates: now, format
using HDF5: HDF5, h5open, create_group, attrs
using Statistics: median
using ShroffCelegansModels
using ShroffCelegansModels: read_config_json, NormalizedDataset,
                             named_cell_orientation_series, resolve_cpaaaa_key,
                             resolve_mapped_key, fix_ap_axis_selectable

# Cells found by survey_cell_orientation_candidates.jl's greedy set-cover
# (min_consistency=0.7, min_n=5) to together cover every strain, using
# strict cell_key.mapping-only resolution. See that script's CSV output for
# the full ranked table this was chosen from.
const DV_CANDIDATE_CELLS = [
    "hyp6_ABplaappap", "DB5", "int7l", "mc2DR", "AVDL", "P12", "G1", "NR_9",
    "DAPPPA_DL11_DBW_Muscle", "URYDR", "ALA_RMED", "hyp7_Caaaaa", "pm3DL", "PDA",
]
const LR_CANDIDATE_CELLS = [
    "P9/10L", "mc2DL", "RIBL", "AVDL", "int1vr", "hyp7_ABarpaappp", "IL1L", "RMDR",
    "DPAAA_VR06_VBW_Muscle", "XXXL", "CANL", "RIPL", "NR_3", "AVKR", "P9/10l", "DD3",
]

named_resolver(target::AbstractString) = (dict, cell_key) -> resolve_mapped_key(dict, cell_key, target)

# (name, axis, resolve) for every check. `hyp7_Cpaaaa` keeps the original,
# more permissive resolver so its coverage isn't narrowed by this change;
# every survey-derived cell uses the strict resolver it was validated with.
const CHECKS = vcat(
    [(name = "hyp7_Cpaaaa", axis = :dv, resolve = resolve_cpaaaa_key)],
    [(name = c, axis = :dv, resolve = named_resolver(c)) for c in DV_CANDIDATE_CELLS],
    [(name = c, axis = :lr, resolve = named_resolver(c)) for c in LR_CANDIDATE_CELLS],
)

axis_sign(a, axis::Symbol) = axis === :dv ? a.dv_sign : a.lr_sign
axis_magnitude(a, axis::Symbol) = axis === :dv ? a.dv_magnitude : a.lr_magnitude

# Majority value among {+1.0, -1.0}, ignoring NaN; NaN if nothing valid.
# Ties default to +1.0 — arbitrary but deterministic.
function majority_sign(xs)
    valid = filter(!isnan, xs)
    isempty(valid) && return NaN
    pos = count(==(1.0), valid)
    neg = count(==(-1.0), valid)
    return pos >= neg ? 1.0 : -1.0
end

# Median of the non-NaN values, or NaN if none.
function robust_median(xs)
    valid = filter(!isnan, xs)
    isempty(valid) && return NaN
    return median(valid)
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

# One check, evaluated across every dataset: per-group per-dataset axes
# series, per-dataset representative sign, and the check's global reference
# sign (majority of representative signs, one dataset one vote).
function evaluate_check(chk, datasets::Dict{String, Vector{NormalizedDataset}})
    series = Dict(
        group => [named_cell_orientation_series(ds, chk.resolve) for ds in dsv]
        for (group, dsv) in datasets
    )
    representative = Dict(
        group => [majority_sign([axis_sign(a, chk.axis) for a in axes]) for axes in series[group]]
        for group in keys(series)
    )
    all_representative = reduce(vcat, values(representative); init=Float64[])
    reference_sign = majority_sign(all_representative)
    return (check = chk, series = series, representative = representative, reference_sign = reference_sign)
end

function save_lattice_orientation(
    filepath::String,
    datasets::Dict{String, Vector{NormalizedDataset}},
    evaluated,
)
    h5open(filepath, "w") do h5f
        A = attrs(h5f)
        A["generated_at"] = format(now(), "yyyy-mm-ddTHH:MM:SS")
        kind_group = create_group(h5f, "lattice_orientation")
        for group in keys(datasets)
            h5g = create_group(kind_group, group)
            for (k, ds) in enumerate(datasets[group])
                dsg = create_group(h5g, string(k))
                dsattrs = attrs(dsg)
                dsattrs["path"] = ds.path
                dsattrs["cell_key.name"] = ds.cell_key.name
                dsattrs["cell_key.start"] = ds.cell_key.start
                dsattrs["cell_key.end"] = ds.cell_key.stop
                dsattrs["cell_key.outliers"] = ds.cell_key.outliers

                check_idx = 0
                for ev in evaluated
                    rep = ev.representative[group][k]
                    isnan(rep) && continue  # this check doesn't resolve for this dataset at all

                    axes = ev.series[group][k]
                    signs = [axis_sign(a, ev.check.axis) for a in axes]
                    magnitudes = [axis_magnitude(a, ev.check.axis) for a in axes]
                    mismatched_timepoints = count(s -> !isnan(s) && s != rep, signs)

                    check_idx += 1
                    cg = create_group(dsg, string("check_", check_idx))
                    cg["sign"] = signs
                    cattrs = attrs(cg)
                    cattrs["annotation_name"] = ev.check.name
                    cattrs["axis"] = string(ev.check.axis)
                    cattrs["magnitude"] = magnitudes
                    cattrs["magnitude_median"] = robust_median(magnitudes)
                    cattrs["representative_sign"] = rep
                    cattrs["reference_sign"] = ev.reference_sign
                    cattrs["matches_reference"] = rep == ev.reference_sign
                    cattrs["mismatched_timepoint_count"] = mismatched_timepoints
                    cattrs["link_eligible"] = fix_ap_axis_selectable(ds.cell_key, ev.check.name)
                end
                dsattrs["check_count"] = check_idx
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

    @info "Checking lattice orientation" check_count=length(CHECKS)
    evaluated = map(chk -> evaluate_check(chk, datasets), CHECKS)
    for ev in evaluated
        all_representative = reduce(vcat, values(ev.representative); init=Float64[])
        resolved = count(!isnan, all_representative)
        mismatched = count(x -> !isnan(x) && x != ev.reference_sign, all_representative)
        @info "Check complete" name=ev.check.name axis=ev.check.axis reference_sign=ev.reference_sign resolved_dataset_count=resolved mismatched_dataset_count=mismatched
    end

    timestamp = format(now(), "yyyy_mm_dd_HHMMSS")
    snapshot_path = joinpath(output_dir, "lattice_orientation_$(timestamp).h5")
    latest_path = joinpath(output_dir, "lattice_orientation.h5")

    @info "Writing snapshot" snapshot_path
    write_atomic(snapshot_path) do tmp
        save_lattice_orientation(tmp, datasets, evaluated)
    end

    @info "Updating latest" latest_path
    write_atomic(latest_path) do tmp
        save_lattice_orientation(tmp, datasets, evaluated)
    end

    @info "Done"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
