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

using Dates: Dates, now, format
using DataFrames: DataFrame
using CSV: CSV
using Printf: @sprintf
using LinearAlgebra: BLAS
using HDF5: h5open, attrs, read_attribute, create_group
using SHA: sha1
using GeometryBasics: Point3

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
    checkpoint_dir = joinpath(output_dir, "checkpoint")
    mkpath(checkpoint_dir)
    ts = format(now(), "yyyy_mm_dd_HHMMSS")
    phase_timings = Pair{String, Float64}[]
    # `body` is first so callers can use do-block syntax: `_phase("...") do ... end`.
    function _phase(body, label::String)
        t = time()
        result = body()
        elapsed_s = round(time() - t; digits=2)
        push!(phase_timings, label => elapsed_s)
        @info "Step done" step=label elapsed_s
        return result
    end

    @info "Pipeline threading" julia_nthreads=Threads.nthreads() blas_nthreads=BLAS.get_num_threads()
    @info "Pipeline starting" config_path output_dir checkpoint_dir n_timepoints kinds

    # 0. Archive any previous run's outputs at the top level into a
    #    timestamped subdir so this run starts with a clean directory.
    #    Top-level always reflects "latest"; consumers (`_latest_cache_path`,
    #    `load_latest_average_annotations`) don't need to know about the archive.
    _phase("0/8 archive_previous") do
        archived = _archive_previous_outputs!(output_dir)
        if archived === nothing
            @info "Nothing to archive" output_dir
        else
            @info "Archived previous outputs" archive=archived
        end
    end

    # 1. Load datasets.
    datasets, flattened = _phase("1/8 load_datasets") do
        @info "[1/8] Loading datasets" config_path
        _, _, datasets = read_config_json(config_path)
        flattened = collect(Iterators.flatten(values(datasets)))
        @info "Loaded datasets" groups=length(datasets) total=length(flattened)
        return datasets, flattened
    end

    # parse_worm_dataset_path.jl pre-loads `annotations_cache` and
    # `my_annotation_position_cache` at module init from HDF5 snapshots baked
    # into the container image — those entries are stale relative to /nearline.
    # Selectively invalidate only datasets whose mtimes have advanced since
    # the cache was populated; unchanged datasets stay cached so the priming
    # step is a no-op for them.
    _phase("1b/8 invalidate_stale") do
        _invalidate_stale_caches!(flattened, kinds)
    end

    # 1c. Load any per-dataset checkpoint files from a previous interrupted
    #     run. Files for paths just invalidated above are skipped.
    n_checkpoints_loaded = _phase("1c/8 load_checkpoints") do
        _load_dataset_checkpoints!(checkpoint_dir, flattened)
    end
    @info "Checkpoints loaded" n=n_checkpoints_loaded dir=checkpoint_dir

    # 2. Prime the annotations_cache for every dataset so subsequent
    #    update_annotations_cache calls have keys to look up. Failures are
    #    logged but non-fatal.
    _phase("2/8 prime_cache") do
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
    end

    # 3 & 4. Load + apply annotation edits.
    n_changes = _phase("3-4/8 load_apply_changes") do
        if isfile(annotation_changes_path)
            @info "[3/8] Loading annotation changes" annotation_changes_path
            changes = ShroffCelegansModels.load_annotation_changes_cache(annotation_changes_path)
            n = length(changes)
            @info "[4/8] Applying annotation changes" n_changes=n
            ShroffCelegansModels.update_annotations_cache(
                ShroffCelegansModels.annotations_cache, changes
            )
            return n
        else
            @warn "Annotation changes file not found — skipping edit application" annotation_changes_path
            return 0
        end
    end

    # 4b. Dated embryo position-cache snapshots (`embryos_*_<date>.h5`,
    #     4 variants: raw/371 × expanded/not). Runs against the
    #     freshly-updated annotations_cache; `clear=false` preserves it.
    #     Failures are non-fatal — the main pipeline outputs still produce.
    _phase("4b/8 save_embryos_dated") do
        @info "[4b/8] Saving dated embryo position cache files" output_dir
        try
            ShroffCelegansModels.save_annotation_position_cache_all_dated(
                datasets;
                clear = false,
                output_dir = output_dir,
            )
        catch err
            @warn "save_annotation_position_cache_all_dated failed (continuing)" err
        end
    end

    # 5. Average lattice models. Disk-cached keyed by (n_timepoints,
    #    max lattice mtime across all datasets) — reuses the existing
    #    save_avg_models / load_avg_models pair, just adds a
    #    `lattice_mtime` HDF5 attribute so we can detect staleness.
    #    Saves ~40s for n=51, more for n=371 when no lattice CSV has
    #    been touched since the cache was written.
    avg_models_h5 = joinpath(output_dir, "avg_models_n$(n_timepoints).h5")
    current_lattice_mtime = _max_lattice_mtime(flattened)
    avg_models = _phase("5/8 get_avg_models") do
        cached = _load_avg_models_if_fresh(avg_models_h5, current_lattice_mtime)
        if cached !== nothing
            @info "[5/8] Loaded cached avg_models" path=avg_models_h5 n_models=length(cached)
            return cached
        end
        @info "[5/8] Computing average lattice models" n_timepoints
        models = get_avg_models(n_timepoints)
        @info "avg_models done" n_models=length(models)
        try
            ShroffCelegansModels.save_avg_models(avg_models_h5, models)
            h5open(avg_models_h5, "r+") do h5f
                attrs(h5f)["lattice_mtime"] = current_lattice_mtime
            end
            @info "Saved avg_models cache" path=avg_models_h5 lattice_mtime=current_lattice_mtime
        catch err
            @warn "avg_models cache save failed (run still valid)" err
        end
        return models
    end

    # 6. Edited averaged annotations.
    avg_dict = _phase("6/8 average_annotations") do
        @info "[6/8] Averaging annotations against avg_models" n_timepoints checkpoint_dir
        ShroffCelegansModels.average_annotations(
            datasets;
            timepoints = LinRange(0, 1, n_timepoints),
            avg_models = avg_models,
            use_cell_key_annotations_only = true,
            checkpoint_dir = checkpoint_dir,
        )
    end

    smoothed = _phase("6.5/8 smooth") do
        @info "[6.5/8] Smoothing" smooth_factor_r smooth_factor_θ smooth_factor_z
        ShroffCelegansModels.smooth_average_annotations(
            avg_dict;
            smooth_factor_r = smooth_factor_r,
            smooth_factor_θ = smooth_factor_θ,
            smooth_factor_z = smooth_factor_z,
        )
    end

    # 7. Write HDF5 in the format the meshscatter web app loads.
    h5_path = joinpath(
        output_dir,
        "edited_smoothed_average_annotations_r$(_factor_token(smooth_factor_r))_theta$(_factor_token(smooth_factor_θ))_z$(_factor_token(smooth_factor_z))_$(ts).h5",
    )
    _phase("7/8 write_h5") do
        @info "[7/8] Writing averaged HDF5" h5_path
        ShroffCelegansModels.save_average_annotations(smoothed; filename = h5_path)
    end

    # 8. Export DataFrame.
    csv_path = joinpath(output_dir, "post_pretwitch_export_$(ts).csv")
    _phase("8/8 export_csv") do
        posttwitch_only = isempty(pretwitch_csv_path) || !isfile(pretwitch_csv_path)
        @info "[8/8] Exporting DataFrame" csv_path pretwitch_csv_path posttwitch_only
        _export_post_pretwitch_csv(
            h5_path, csv_path;
            pretwitch_csv_path = posttwitch_only ? nothing : pretwitch_csv_path,
            avg_models = avg_models,
        )
    end

    # 8c. Combined pretwitch + posttwitch explicit-schema CSVs. Uses the
    #     `_for_ben.csv` produced by `resave_for_ben` (called from phase 8)
    #     plus the ryan_data pretwitch coords + naming-correlations table
    #     baked into the repo. Positional model cell names are translated
    #     to embryonic lineage names.
    ben_csv = replace(h5_path, ".h5" => "_for_ben.csv")
    _phase("8a/8 export_combined") do
        if isfile(ben_csv)
            try
                ShroffCelegansModels.write_combined_explicit_csvs(;
                    output_dir = output_dir,
                    avg_models = avg_models,
                    ben_csv_path = ben_csv,
                    date_str = ts,
                )
                @info "[8c/8] Wrote combined pretwitch+posttwitch CSVs"
            catch err
                @warn "Combined explicit-CSV export failed (pipeline outputs still valid)" err
            end
        else
            @warn "Skipping 8c — _for_ben.csv not found" ben_csv
        end
    end

    # Persist the in-memory caches so future package boots (interactive
    # sessions, web service restarts) load the freshly-recomputed data
    # instead of the stale snapshots baked into the container image.
    # parse_worm_dataset_path.jl picks up these files at module init when
    # `RECOMPUTE_OUTPUT_DIR` env is set or its default exists.
    annotations_cache_h5 = joinpath(output_dir, "annotations_cache.h5")
    my_positions_h5 = joinpath(output_dir, "my_annotation_position_cache.h5")
    _phase("8b/8 persist_caches") do
        try
            _write_atomic(annotations_cache_h5) do tmp
                ShroffCelegansModels.save_annotations_cache(ShroffCelegansModels.annotations_cache; filename = tmp)
            end
            _write_atomic(my_positions_h5) do tmp
                ShroffCelegansModels.save_annotation_cache(; filename = tmp)
            end
            @info "Persisted caches" annotations_cache_h5 my_positions_h5
        catch err
            @warn "Cache persistence failed (pipeline outputs still valid)" err
        end
    end

    # Write a static index.html so the /recompute/ route shows full
    # filenames + sizes + mtimes (nginx's autoindex truncates at ~50
    # chars). Generates one for the top-level dir and one for each
    # `archive_*/` subdir.
    _phase("8c/8 write_index") do
        try
            _write_recompute_index(output_dir)
            @info "Wrote recompute index" output_dir
        catch err
            @warn "Index generation failed (pipeline outputs still valid)" err
        end
    end

    # On full success, delete the checkpoint dir so the next run starts fresh.
    try
        rm(checkpoint_dir; recursive=true, force=true)
    catch err
        @warn "Could not remove checkpoint dir after success" checkpoint_dir err
    end

    @info "Pipeline complete" h5_path csv_path n_changes n_timepoints phase_timings
    return (; h5_path, csv_path, n_changes, n_timepoints, phase_timings, annotations_cache_h5, my_positions_h5)
end

# Sweep dated top-level outputs into `archive_<ts>/`, where `<ts>` is the
# newest run-timestamp embedded in any existing filename. Used for
# per-run products (edited_smoothed_*_<ts>.h5, post_pretwitch_export_*.csv,
# embryos_*_<date>.h5, etc.).
#
# Left in place:
#   - subdirectories (already-archived, `checkpoint/`, etc.)
#   - in-progress `.tmp.*` writes
#   - "rolling" cache files that get overwritten each run (no date in name):
#     annotations_cache.h5, my_annotation_position_cache.h5, avg_models_n*.h5
#
# A file is considered "dated" if it contains a `YYYY_MM_DD` substring (with
# or without the trailing `_HHMMSS`). Cache files have no date and therefore
# stay at the top level for the next run to reuse.
#
# Returns the archive dir path, or `nothing` if there's nothing to archive.
const _DATE_RE = r"\d{4}_\d{2}_\d{2}"
const _TS_RE   = r"(\d{4}_\d{2}_\d{2}_\d{6})"

function _archive_previous_outputs!(output_dir::AbstractString)
    isdir(output_dir) || return nothing

    latest_ts = ""
    candidates = String[]
    for f in readdir(output_dir)
        full = joinpath(output_dir, f)
        isdir(full) && continue              # skip checkpoint/, archive_*/, ...
        occursin(".tmp.", f) && continue     # skip in-progress writes
        occursin(_DATE_RE, f) || continue    # skip rolling caches (no date)

        push!(candidates, f)
        m = match(_TS_RE, f)
        if m !== nothing && m.captures[1] > latest_ts
            latest_ts = m.captures[1]
        end
    end

    isempty(candidates) && return nothing
    if isempty(latest_ts)
        # Only dated-but-not-full-timestamp files (e.g. only embryos_*_<date>.h5);
        # fall back to the date portion of the first candidate.
        m = match(_DATE_RE, first(candidates))
        latest_ts = m === nothing ? "previous" : m.match
    end

    archive_dir = joinpath(output_dir, "archive_" * latest_ts)
    if isdir(archive_dir)
        @warn "Archive dir already exists; skipping move to avoid overwrite" archive_dir
        return nothing
    end

    mkpath(archive_dir)
    n_moved = 0
    for f in candidates
        src = joinpath(output_dir, f)
        dst = joinpath(archive_dir, f)
        try
            mv(src, dst)
            n_moved += 1
        catch err
            @warn "Failed to archive file" src dst err
        end
    end
    @info "Archive populated" archive_dir n_moved of=length(candidates)
    return archive_dir
end

# Max lattice mtime across every dataset's every timepoint — the staleness
# signal for the avg_models disk cache. avg_models depends only on lattice
# data (it averages CelegansModels), so a single global max across all
# datasets is sufficient: if any lattice CSV is newer than the cache, we
# recompute.
function _max_lattice_mtime(datasets::Vector{<:ShroffCelegansModels.Datasets.NormalizedDataset})::Float64
    m = NaN
    for ds in datasets
        for x in ShroffCelegansModels.MIPAVIO.get_lattice_modified_times_unix(ds)
            isnan(x) && continue
            if isnan(m) || x > m
                m = x
            end
        end
    end
    return m
end

# Load an avg_models HDF5 written by `save_avg_models` IF it carries a
# `lattice_mtime` attr that's at least `current_lattice_mtime`. Returns
# `nothing` for missing file, missing attr, or stale cache.
function _load_avg_models_if_fresh(path::AbstractString, current_lattice_mtime::Float64)
    isfile(path) || return nothing
    isnan(current_lattice_mtime) && return nothing
    stored_mtime = try
        h5open(path, "r") do h5f
            if haskey(attrs(h5f), "lattice_mtime")
                Float64(read(attrs(h5f)["lattice_mtime"]))
            else
                NaN
            end
        end
    catch err
        @warn "Could not read lattice_mtime from avg_models cache" path err
        return nothing
    end
    if isnan(stored_mtime) || stored_mtime < current_lattice_mtime
        @info "avg_models cache stale" path stored_mtime current_lattice_mtime
        return nothing
    end
    return try
        ShroffCelegansModels.load_avg_models(path)
    catch err
        @warn "Failed to load avg_models cache; will recompute" path err
        nothing
    end
end

# Atomic write: invoke `body(tmp_path)` to produce the file, then mv it into
# place. POSIX rename is atomic on the same filesystem; readers always see
# either the previous version or the new one — never a partial write.
function _write_atomic(body, dst::AbstractString)
    tmp = string(dst, ".tmp.", getpid(), ".", time_ns())
    try
        body(tmp)
        mv(tmp, dst; force=true)
    catch
        isfile(tmp) && rm(tmp; force=true)
        rethrow()
    end
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

# Convert a Windows-style cache key path ("X:\foo\bar") to the equivalent
# Linux path under /nearline/shroff. Legacy `annotations_cache.h5` entries
# baked into the container image use Windows separators; the live cache
# (populated by load_straightened_annotations_over_time) uses Linux paths
# matching `dataset.path`. Normalizing the cache key lets us match either
# format against `dataset.path` for staleness comparison.
function _normalize_cache_path(p::AbstractString)::String
    s = String(p)
    occursin('\\', s) || return s
    # Drive-letter prefix `X:\foo\bar` → `/nearline/shroff/foo/bar`.
    if length(s) >= 3 && isuppercase(s[1]) && s[2] == ':' && s[3] == '\\'
        s = "/nearline/shroff/" * s[4:end]
    end
    return replace(s, "\\" => "/")
end

# Selective cache invalidation. Stats current annotation+lattice mtimes for
# every dataset, then removes any `annotations_cache` entries whose cached
# mtime predates the current disk state (or is NaN — i.e. legacy entries
# loaded from the pre-mtime-tracking HDF5 schema).
#
# For `my_annotation_position_cache` (which stores positions transformed via
# avg_models and has no mtime tracking): if `kinds` contains "lattice", clear
# all entries because avg_models change globally; otherwise drop just the
# entries whose corresponding annotations_cache key was invalidated.
function _invalidate_stale_caches!(
    datasets::Vector{<:ShroffCelegansModels.Datasets.NormalizedDataset},
    kinds::Vector{String},
)
    # Per-dataset max mtime across annotation + lattice inputs. Mirrors what
    # `_dataset_mtime` records into `AnnotationsCacheValue.mtime`.
    current_max = Dict{String, Float64}()
    for ds in datasets
        ann = ShroffCelegansModels.MIPAVIO.get_annotation_modified_times_unix(ds)
        lat = ShroffCelegansModels.MIPAVIO.get_lattice_modified_times_unix(ds)
        all_m = Float64[]
        for m in ann; isnan(m) || push!(all_m, m); end
        for m in lat; isnan(m) || push!(all_m, m); end
        current_max[ds.path] = isempty(all_m) ? NaN : maximum(all_m)
    end

    annotations_cache = ShroffCelegansModels.annotations_cache
    n_before = length(annotations_cache)
    invalidated_paths = Set{String}()
    unmatched_paths = Set{String}()
    for k in collect(keys(annotations_cache))
        cached_path = k[1]
        norm = _normalize_cache_path(cached_path)
        if !haskey(current_max, norm)
            # Cache entry references a dataset that's not in the current config.
            # Conservative: leave it alone (might be referenced by an interactive
            # session). Track for visibility.
            push!(unmatched_paths, norm)
            continue
        end
        cur = current_max[norm]
        cached_mt = annotations_cache[k].mtime
        is_stale = isnan(cached_mt) || (!isnan(cur) && cur > cached_mt)
        if is_stale
            delete!(annotations_cache, k)
            push!(invalidated_paths, norm)
        end
    end
    @info "annotations_cache invalidation" n_before invalidated=length(invalidated_paths) remaining=length(annotations_cache) unmatched=length(unmatched_paths)

    my_cache = ShroffCelegansModels.my_annotation_position_cache
    if "lattice" in kinds
        n_my_before = length(my_cache)
        empty!(my_cache)
        @info "my_annotation_position_cache fully cleared (lattice change → avg_models will be recomputed)" n_my_before
    else
        n_removed = 0
        for p in invalidated_paths
            if haskey(my_cache, p)
                delete!(my_cache, p)
                n_removed += 1
            end
        end
        @info "my_annotation_position_cache selective invalidation" n_removed remaining=length(my_cache)
    end
end

# 16-hex-char filename derived from sha1(path) — enough collision resistance
# for our ~72-dataset universe and short enough to fit comfortably in any FS.
checkpoint_filename(path::AbstractString) = bytes2hex(sha1(path))[1:16] * ".h5"

# Write a single dataset's positions to `<checkpoint_dir>/<hash>.h5` atomically
# (via temp + rename). Stores the full `dataset.path` as an attribute for
# inspection. Returns the bytes written (0 on no-op). Called from the threaded
# loop in get_group_annotation_positions_over_time; no lock needed because each
# dataset has a unique filename.
function write_dataset_checkpoint(
    checkpoint_dir::AbstractString,
    dataset_path::AbstractString,
    positions::Vector{Vector{Point3{Float64}}},
)::Int
    isempty(checkpoint_dir) && return 0
    final_path = joinpath(checkpoint_dir, checkpoint_filename(dataset_path))
    tmp_path = string(final_path, ".tmp.", getpid(), ".", time_ns())
    try
        h5open(tmp_path, "w") do h5f
            attrs(h5f)["path"] = String(dataset_path)
            g = create_group(h5f, "positions")
            for (idx, pts) in pairs(positions)
                # Layout matches save_annotation_cache: 3×N row-stacked matrix.
                m = reinterpret(Float64, pts)
                m = reshape(m, 3, :)
                tp_name = @sprintf("timepoint_%03d", idx)
                g[tp_name] = collect(transpose(m))
            end
        end
        mv(tmp_path, final_path; force=true)
        return filesize(final_path)
    catch err
        isfile(tmp_path) && rm(tmp_path; force=true)
        @warn "Checkpoint write failed" final_path dataset_path err
        return 0
    end
end

# Read all `<checkpoint_dir>/*.h5` files (skipping stale `*.tmp.*` artefacts)
# and insert each into `my_annotation_position_cache`. Entries whose `path`
# attribute is NOT in the current dataset list are ignored (treated as stale,
# left on disk for human inspection — likely from an older config). Returns
# the count loaded.
function _load_dataset_checkpoints!(
    checkpoint_dir::AbstractString,
    datasets::Vector{<:ShroffCelegansModels.Datasets.NormalizedDataset},
)::Int
    isdir(checkpoint_dir) || return 0
    # Build the valid-paths set from CURRENT datasets, not what was on disk
    # last time — drops cross-config stale checkpoints.
    valid_paths = Set(ds.path for ds in datasets)
    my_cache = ShroffCelegansModels.my_annotation_position_cache
    n_loaded = 0
    n_skipped_invalidated = 0
    n_skipped_other = 0
    for fname in readdir(checkpoint_dir)
        endswith(fname, ".h5") || continue
        occursin(".tmp.", fname) && (rm(joinpath(checkpoint_dir, fname); force=true); continue)
        fpath = joinpath(checkpoint_dir, fname)
        try
            h5open(fpath, "r") do h5f
                ds_path = String(read_attribute(h5f, "path"))
                if !(ds_path in valid_paths)
                    n_skipped_other += 1
                    return
                end
                # If the in-memory `my_annotation_position_cache` was just
                # invalidated for this path (selective invalidation removed it),
                # the checkpoint reinstates it — i.e. resume from where the
                # previous run left off. If it wasn't invalidated, we'd just
                # be loading data identical to what's already cached; still
                # cheap and tolerant.
                g = h5f["positions"]
                tps = sort(parse.(Int, last.(split.(filter(startswith("timepoint_"), keys(g)), "_"))))
                positions = Vector{Vector{Point3{Float64}}}(undef, length(tps))
                for tp_i in tps
                    tp_name = @sprintf("timepoint_%03d", tp_i)
                    mat = transpose(g[tp_name][])::AbstractMatrix{Float64}
                    positions[tp_i] = vec(reinterpret(Point3{Float64}, mat))
                end
                my_cache[ds_path] = positions
                n_loaded += 1
            end
        catch err
            @warn "Failed to load checkpoint file" fpath err
        end
    end
    n_skipped_other > 0 && @info "Checkpoint files skipped (path not in current config)" n=n_skipped_other
    return n_loaded
end

# ---------- static index.html generation for /recompute/ ----------
#
# nginx's autoindex (html mode) truncates filenames to ~50 chars. We
# generate our own index.html that nginx serves automatically when
# present. Pico.css matches the rest of the Shroff web pages.

const _PICO_CSS_LINK = "<link rel=\"stylesheet\" href=\"https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css\">"
const _SHROFF_CSS_LINK = "<link rel=\"stylesheet\" href=\"/style.css\">"

function _human_size(n::Integer)
    n < 1024 && return string(n, " B")
    n < 1024^2 && return string(round(n / 1024; digits=1), " KB")
    n < 1024^3 && return string(round(n / 1024^2; digits=1), " MB")
    return string(round(n / 1024^3; digits=2), " GB")
end

_format_mtime(m::Float64) = isnan(m) ? "" : Dates.format(Dates.unix2datetime(m), "yyyy-mm-dd HH:MM:SS")

# Sanitize a string for use as an HTML text node. We control the inputs
# (filenames, sizes, timestamps) but escape anyway for safety.
function _html_escape(s::AbstractString)
    out = IOBuffer()
    for c in s
        c == '&' ? print(out, "&amp;") :
        c == '<' ? print(out, "&lt;") :
        c == '>' ? print(out, "&gt;") :
        c == '"' ? print(out, "&quot;") :
        c == '\'' ? print(out, "&#39;") :
        print(out, c)
    end
    return String(take!(out))
end

function _render_index_html(;
    title::AbstractString,
    intro::AbstractString,
    files::Vector,                       # Vector of NamedTuples (name, size, mtime)
    subdirs::Vector{<:AbstractString},
    parent_link::AbstractString,
)
    rows = isempty(files) ? "<tr><td colspan=\"3\"><em>(no files)</em></td></tr>" :
        join(map(files) do f
            name_esc = _html_escape(f.name)
            "<tr><td><a href=\"$(name_esc)\"><code>$(name_esc)</code></a></td>" *
            "<td>$(_human_size(f.size))</td>" *
            "<td>$(_format_mtime(f.mtime))</td></tr>"
        end, "\n            ")
    subdir_section = if isempty(subdirs)
        ""
    else
        items = join(map(subdirs) do d
            d_esc = _html_escape(d)
            "<li><a href=\"$(d_esc)/\"><code>$(d_esc)/</code></a></li>"
        end, "\n              ")
        """

            <h2>Archived runs</h2>
            <ul>
              $(items)
            </ul>"""
    end
    generated_at = _format_mtime(Float64(time()))
    title_esc = _html_escape(title)
    intro_esc = _html_escape(intro)
    """<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>$(title_esc)</title>
    $(_PICO_CSS_LINK)
    $(_SHROFF_CSS_LINK)
  </head>
  <body>
    <main class="container">
      <p><a href="$(parent_link)">&larr; up</a></p>
      <h1>$(title_esc)</h1>
      <p>$(intro_esc)</p>
      <p><small>Generated $(generated_at)</small></p>
      <h2>Files</h2>
      <table>
        <thead><tr><th>Name</th><th>Size</th><th>Modified</th></tr></thead>
        <tbody>
            $(rows)
        </tbody>
      </table>$(subdir_section)
    </main>
  </body>
</html>
"""
end

# Walks `output_dir`, writes index.html at top level + one per `archive_*/`.
# Treats `checkpoint/` and any other non-archive subdirs as opaque (no index).
function _write_recompute_index(output_dir::AbstractString)
    isdir(output_dir) || return nothing

    files = String[]
    archive_dirs = String[]
    for entry in readdir(output_dir)
        full = joinpath(output_dir, entry)
        if isdir(full)
            startswith(entry, "archive_") && push!(archive_dirs, entry)
        else
            entry == "index.html" && continue
            occursin(".tmp.", entry) && continue
            push!(files, entry)
        end
    end
    sort!(files, by = f -> mtime(joinpath(output_dir, f)), rev = true)
    sort!(archive_dirs, rev = true)

    top_html = _render_index_html(;
        title = "Recompute outputs",
        intro = "Top-level files are the latest run's outputs (latest cache files have no date in their name). " *
                "Archive subdirectories preserve previous runs.",
        files = [(name = f, size = filesize(joinpath(output_dir, f)), mtime = mtime(joinpath(output_dir, f))) for f in files],
        subdirs = archive_dirs,
        parent_link = "/",
    )
    write(joinpath(output_dir, "index.html"), top_html)

    for ad in archive_dirs
        adir = joinpath(output_dir, ad)
        afiles = String[]
        for entry in readdir(adir)
            entry == "index.html" && continue
            isfile(joinpath(adir, entry)) || continue
            push!(afiles, entry)
        end
        sort!(afiles, by = f -> mtime(joinpath(adir, f)), rev = true)
        archive_html = _render_index_html(;
            title = "Archive: " * ad,
            intro = "Archived recompute outputs from a previous run.",
            files = [(name = f, size = filesize(joinpath(adir, f)), mtime = mtime(joinpath(adir, f))) for f in afiles],
            subdirs = String[],
            parent_link = "../",
        )
        write(joinpath(adir, "index.html"), archive_html)
    end
    return nothing
end
