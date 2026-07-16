"""
Surveys every annotated cell across every dataset in the active config_path
and ranks candidates for use as alternate LR/DV orientation anchors, for
strains where Cpaaaa (hyp7) isn't tracked. Reuses the same untwisted-frame
signal as `check_lattice_orientation.jl`'s Cpaaaa check
(`ShroffCelegansModels.dataset_cell_orientation_survey`), but computed for
every cell name instead of just Cpaaaa/hyp7.

For each cell, across every timepoint of every dataset where it's annotated,
records the DV (dorsal/ventral, `y`) and LR (`x`) sign+magnitude relative to
the seam-cell plane. A good anchor candidate is a cell that is:
- consistently one sign (DV or LR) across the vast majority of timepoints
  it's seen in — `*_consistency` close to 1.0.
- far from the plane on that axis (large `*_magnitude_median`) — a small
  magnitude means the sign flips easily from noise even when technically
  "consistent".
- present in many datasets/strains — a candidate seen in only one or two
  datasets isn't a reliable population-wide anchor.

`*_score = *_consistency * *_magnitude_median` combines the first two into
a single ranking number; presence (`dataset_count`/`group_count`) is
reported alongside so low-presence cells can be filtered by eye even if a
hard `MIN_DATASETS`/`MIN_GROUPS` cutoff (env vars, defaults 5/2) already
excludes the least-supported ones from the ranked printout.

Writes a full CSV (`CELL_ORIENTATION_CANDIDATES_PATH`, default
`cell_orientation_candidates.csv` in the working directory) with one row per
cell name (unfiltered by presence), and logs the top candidates by DV score
and by LR score. This is a one-off/periodic research tool, not part of the
production QC pipeline — run it manually (locally, or via an ad hoc cluster
Job) rather than on a schedule.
"""

using Statistics: median
using ShroffCelegansModels
using ShroffCelegansModels: read_config_json, dataset_cell_orientation_survey

mutable struct CellAgg
    dv_signs::Vector{Float64}
    dv_magnitudes::Vector{Float64}
    lr_signs::Vector{Float64}
    lr_magnitudes::Vector{Float64}
    datasets::Set{String}
    groups::Set{String}
end
CellAgg() = CellAgg(Float64[], Float64[], Float64[], Float64[], Set{String}(), Set{String}())

# Fraction of a sign vector's non-NaN entries taken by the majority sign,
# plus the number of non-NaN entries it was computed from.
function consistency(signs::AbstractVector{<:Real})
    valid = filter(!isnan, signs)
    isempty(valid) && return NaN, 0
    pos = count(==(1.0), valid)
    neg = count(==(-1.0), valid)
    return max(pos, neg) / length(valid), length(valid)
end

robust_median(xs) = (valid = filter(!isnan, xs); isempty(valid) ? NaN : median(valid))

function survey_all(datasets::AbstractDict)
    agg = Dict{String, CellAgg}()
    for (group, group_datasets) in datasets
        for (idx, ds) in enumerate(group_datasets)
            @info "Surveying dataset" group idx path=ds.path
            per_cell = dataset_cell_orientation_survey(ds)
            dataset_id = string(group, "#", idx)
            for (name, axes_vec) in per_cell
                a = get!(agg, name, CellAgg())
                for ax in axes_vec
                    push!(a.dv_signs, ax.dv_sign)
                    push!(a.dv_magnitudes, ax.dv_magnitude)
                    push!(a.lr_signs, ax.lr_sign)
                    push!(a.lr_magnitudes, ax.lr_magnitude)
                end
                push!(a.datasets, dataset_id)
                push!(a.groups, group)
            end
        end
    end
    return agg
end

function summarize(agg::Dict{String, CellAgg})
    rows = map(collect(agg)) do (name, a)
        dv_cons, dv_n = consistency(a.dv_signs)
        lr_cons, lr_n = consistency(a.lr_signs)
        dv_mag = robust_median(a.dv_magnitudes)
        lr_mag = robust_median(a.lr_magnitudes)
        dv_score = (isnan(dv_cons) || isnan(dv_mag)) ? -Inf : dv_cons * dv_mag
        lr_score = (isnan(lr_cons) || isnan(lr_mag)) ? -Inf : lr_cons * lr_mag
        (
            name = name,
            dataset_count = length(a.datasets),
            group_count = length(a.groups),
            dv_consistency = dv_cons,
            dv_n = dv_n,
            dv_magnitude_median = dv_mag,
            dv_score = dv_score,
            lr_consistency = lr_cons,
            lr_n = lr_n,
            lr_magnitude_median = lr_mag,
            lr_score = lr_score,
        )
    end
    return rows
end

function write_csv(path::AbstractString, rows)
    header = [
        "name", "dataset_count", "group_count",
        "dv_consistency", "dv_n", "dv_magnitude_median", "dv_score",
        "lr_consistency", "lr_n", "lr_magnitude_median", "lr_score",
    ]
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join((r.name, r.dataset_count, r.group_count,
                               r.dv_consistency, r.dv_n, r.dv_magnitude_median, r.dv_score,
                               r.lr_consistency, r.lr_n, r.lr_magnitude_median, r.lr_score), ","))
        end
    end
end

function log_top(rows, score_key::Symbol, label::AbstractString; top_n::Int=20)
    ranked = sort(rows; by=r -> getproperty(r, score_key), rev=true)
    @info "Top $label candidates" [
        (name=r.name, score=round(getproperty(r, score_key); digits=3),
         datasets=r.dataset_count, groups=r.group_count)
        for r in first(ranked, min(top_n, length(ranked)))
    ]
end

function main()
    min_datasets = parse(Int, get(ENV, "MIN_DATASETS", "5"))
    min_groups = parse(Int, get(ENV, "MIN_GROUPS", "2"))
    output_path = get(ENV, "CELL_ORIENTATION_CANDIDATES_PATH", "cell_orientation_candidates.csv")

    config_path = ShroffCelegansModels.config_path
    @info "Loading datasets from config" config_path
    _, _, datasets = read_config_json(config_path)
    @info "Loaded datasets" groups=length(datasets) total=sum(length, values(datasets))

    agg = survey_all(datasets)
    rows = summarize(agg)

    @info "Writing candidate table" output_path cell_count=length(rows)
    write_csv(output_path, rows)

    supported = filter(r -> r.dataset_count >= min_datasets && r.group_count >= min_groups, rows)
    @info "Cells meeting presence threshold" min_datasets min_groups count=length(supported)
    log_top(supported, :dv_score, "DV (dorsal/ventral)")
    log_top(supported, :lr_score, "LR (left/right)")

    @info "Done"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
