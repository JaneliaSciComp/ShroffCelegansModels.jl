"""
End-to-end recompute pipeline. Replaces the stub `run_pipeline` previously
embedded in `web/scripts/run_recompute_if_needed.jl`.

Orchestrates the user's 6-step workflow:

  1. Load datasets from the active config.
  2. Prime the in-memory `annotations_cache` so subsequent edits can be applied
     against existing keys (per the workflow note: changes can't apply against
     a cold cache).
  3. Load annotation edits from `annotation_changes.h5` (the file maintained by
     the interactive web service).
  4. Apply the edits via `update_annotations_cache`.
  5. Compute average lattice models via `get_avg_models(n_timepoints)`.
  6. Compute edited averaged annotations via `average_annotations(...)`, then
     smooth via `smooth_average_annotations(...)`.
  7. Write the smoothed averages to a timestamped HDF5 file.
  8. Export the post-twitch DataFrame (and pre-twitch if the source CSV is
     available) to CSV.

The function takes no required arguments — it reads its configuration from env
vars / package defaults — so an LSF wrapper could invoke it identically once
/nearline is bridged elsewhere. No OpenShift-specific assumptions live here.
"""

using Dates: now, format
using DataFrames: DataFrame
using CSV: CSV
using Printf: @sprintf

# Default smoothing matches the production filename pattern
# `..._r020_theta020_z030_...`.
const _DEFAULT_SMOOTH_R = 0.20
const _DEFAULT_SMOOTH_θ = 0.20
const _DEFAULT_SMOOTH_Z = 0.30

"""
    run_recompute_pipeline(; kwargs...) -> NamedTuple

Run the full recompute pipeline. Returns a NamedTuple of artifact paths and
counts.

Keyword arguments (all with sensible defaults):

  - `config_path` — path to the datasets config JSON
  - `output_dir` — where to write artifacts (HDF5 + CSV). Created if absent.
  - `annotation_changes_path` — input edits HDF5. Skip step 3/4 if missing.
  - `pretwitch_csv_path` — optional pre-twitch coordinates CSV. If absent the
    DataFrame export contains post-twitch data only.
  - `n_timepoints` — `LinRange(0, 1, n_timepoints)` for averaging. Default 371.
  - `kinds` — sorted-unique kinds reported in the trigger marker. Currently
    informational; future versions may use it to skip get_avg_models when only
    annotations changed.
  - `smooth_factor_r`, `smooth_factor_θ`, `smooth_factor_z` — passed to
    `smooth_average_annotations`. Defaults match the production output
    `..._r020_theta020_z030_...`.
"""
function run_recompute_pipeline(;
    config_path::AbstractString = ShroffCelegansModels.config_path,
    output_dir::AbstractString = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute"),
    annotation_changes_path::AbstractString = get(ENV, "ANNOTATION_CHANGES_PATH", "/data/annotations/annotation_changes.h5"),
    pretwitch_csv_path::AbstractString = get(ENV, "PRETWITCH_CSV_PATH", ""),
    n_timepoints::Int = 371,
    kinds::Vector{String} = ["annotation", "lattice"],
    smooth_factor_r::Float64 = _DEFAULT_SMOOTH_R,
    smooth_factor_θ::Float64 = _DEFAULT_SMOOTH_θ,
    smooth_factor_z::Float64 = _DEFAULT_SMOOTH_Z,
)
    mkpath(output_dir)
    ts = format(now(), "yyyy_mm_dd_HHMMSS")
    @info "Pipeline starting" config_path output_dir n_timepoints kinds

    # Clear pre-loaded caches before doing anything else. parse_worm_dataset_path.jl
    # populates `annotations_cache` and `my_annotation_position_cache` at module
    # init from baked-in HDF5 snapshots (annotations_cache.h5 +
    # my_annotation_position_cache.h5). Those snapshots are stale relative to
    # /nearline — the whole point of this pipeline is to recompute against
    # current disk state, so we must blow away the pre-loaded entries first.
    # Without this, the priming step and average_annotations would short-circuit
    # on `haskey(cache, ...)` and return stale data.
    @info "Clearing pre-loaded caches before recompute"
    empty!(ShroffCelegansModels.annotations_cache)
    empty!(ShroffCelegansModels.my_annotation_position_cache)

    # 1. Load datasets.
    @info "[1/8] Loading datasets" config_path
    _, _, datasets = read_config_json(config_path)
    flattened = collect(Iterators.flatten(values(datasets)))
    @info "Loaded datasets" groups=length(datasets) total=length(flattened)

    # 2. Prime the annotations_cache for every dataset so subsequent
    #    update_annotations_cache calls have keys to look up. Failures are
    #    logged but non-fatal.
    @info "[2/8] Priming annotations_cache"
    primed = 0
    for ds in flattened
        try
            load_straightened_annotations_over_time(ds; use_myuntwist=true)
            primed += 1
        catch err
            @warn "Cache prime failed for dataset" path=ds.path err
        end
    end
    @info "Cache primed" datasets=primed of=length(flattened)

    # 3 & 4. Load + apply annotation edits.
    n_changes = 0
    if isfile(annotation_changes_path)
        @info "[3/8] Loading annotation changes" annotation_changes_path
        changes = ShroffCelegansModels.load_annotation_changes_cache(annotation_changes_path)
        n_changes = length(changes)
        @info "[4/8] Applying annotation changes" n_changes
        ShroffCelegansModels.update_annotations_cache(
            ShroffCelegansModels.annotations_cache, changes
        )
    else
        @warn "Annotation changes file not found — skipping edit application" annotation_changes_path
    end

    # 5. Average lattice models.
    @info "[5/8] Computing average lattice models" n_timepoints
    t0 = time()
    avg_models = get_avg_models(n_timepoints)
    @info "avg_models done" elapsed_s=round(time() - t0; digits=1) n_models=length(avg_models)

    # 6. Edited averaged annotations + smoothing.
    @info "[6/8] Averaging annotations against avg_models" n_timepoints
    t0 = time()
    avg_dict = ShroffCelegansModels.average_annotations(
        datasets;
        timepoints = LinRange(0, 1, n_timepoints),
        avg_models = avg_models,
        use_cell_key_annotations_only = true,
    )
    @info "average_annotations done" elapsed_s=round(time() - t0; digits=1)

    @info "[6.5/8] Smoothing" smooth_factor_r smooth_factor_θ smooth_factor_z
    smoothed = ShroffCelegansModels.smooth_average_annotations(
        avg_dict;
        smooth_factor_r = smooth_factor_r,
        smooth_factor_θ = smooth_factor_θ,
        smooth_factor_z = smooth_factor_z,
    )

    # 7. Write HDF5 in the format the meshscatter web app loads.
    h5_path = joinpath(
        output_dir,
        "edited_smoothed_average_annotations_r$(_factor_token(smooth_factor_r))_theta$(_factor_token(smooth_factor_θ))_z$(_factor_token(smooth_factor_z))_$(ts).h5",
    )
    @info "[7/8] Writing averaged HDF5" h5_path
    ShroffCelegansModels.save_average_annotations(smoothed; filename = h5_path)

    # 8. Export DataFrame. Always produce the post-twitch CSV; combine with
    #    pre-twitch if the source file is configured and present.
    csv_path = joinpath(output_dir, "post_pretwitch_export_$(ts).csv")
    posttwitch_only = isempty(pretwitch_csv_path) || !isfile(pretwitch_csv_path)
    @info "[8/8] Exporting DataFrame" csv_path pretwitch_csv_path posttwitch_only
    _export_post_pretwitch_csv(
        h5_path, csv_path;
        pretwitch_csv_path = posttwitch_only ? nothing : pretwitch_csv_path,
        avg_models = avg_models,
    )

    @info "Pipeline complete" h5_path csv_path n_changes n_timepoints
    return (; h5_path, csv_path, n_changes, n_timepoints)
end

# "0.20" -> "020"; "0.3" -> "030". Matches the existing on-disk filename
# convention `..._r020_theta020_z030_...`.
function _factor_token(x::Float64)
    s = @sprintf("%03d", round(Int, x * 100))
    return s
end

# Wrap resave_for_ben to produce the post-twitch CSV, optionally vcat'ing the
# pre-twitch source CSV (translated to the same column schema). If no pretwitch
# source is provided/found, the output is post-twitch only — non-fatal.
function _export_post_pretwitch_csv(
    h5_path::AbstractString,
    output_csv::AbstractString;
    pretwitch_csv_path::Union{Nothing, AbstractString},
    avg_models,
)
    # resave_for_ben writes <basename>_for_ben.csv plus _ryan_duplicates.csv
    # and _ryan_stats.csv. The Ben CSV has columns (cell, time, x, y, z) over
    # the post-twitch window 381–751 mpfc.
    ben_csv = replace(h5_path, ".h5" => "_for_ben.csv")
    isfile(ben_csv) && rm(ben_csv)
    ben_aux1 = replace(ben_csv, ".csv" => "_ryan_duplicates.csv")
    isfile(ben_aux1) && rm(ben_aux1)
    ben_aux2 = replace(ben_csv, ".csv" => "_ryan_stats.csv")
    isfile(ben_aux2) && rm(ben_aux2)

    ShroffCelegansModels.resave_for_ben(h5_path; target_filename = ben_csv, time_range = (381, 751))
    posttwitch_df = CSV.read(ben_csv, DataFrame)

    # Normalize to the explicit-export schema:
    # (lineage_name, minutes_post_first_cleavage, LR_micrometers, DV_micrometers, AP_micrometers)
    posttwitch_explicit = DataFrame(
        lineage_name = posttwitch_df.cell,
        minutes_post_first_cleavage = posttwitch_df.time,
        LR_micrometers = posttwitch_df.x,
        DV_micrometers = posttwitch_df.z,
        AP_micrometers = posttwitch_df.y,
    )

    if pretwitch_csv_path === nothing
        CSV.write(output_csv, posttwitch_explicit; writeheader = true)
        return output_csv
    end

    # Pre-twitch source: pre-existing user file with at minimum
    # (cell, time, x, y, z) columns matching the historical format. We trust
    # it and just rename columns into the explicit schema; downstream consumers
    # can post-process if needed.
    pretwitch_df = CSV.read(pretwitch_csv_path, DataFrame)
    pretwitch_explicit = DataFrame(
        lineage_name = pretwitch_df.cell,
        minutes_post_first_cleavage = pretwitch_df.time,
        LR_micrometers = pretwitch_df.z,
        DV_micrometers = pretwitch_df.y,
        AP_micrometers = pretwitch_df.x,
    )

    combined = vcat(pretwitch_explicit, posttwitch_explicit)
    CSV.write(output_csv, combined; writeheader = true)
    return output_csv
end
