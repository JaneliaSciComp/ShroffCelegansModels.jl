"""
Runs the recompute pipeline if the daily `save_modified_times.jl` CronJob
detected mtime changes since its previous snapshot.

Pattern: filesystem marker on the shared annotations PVC.

  - `MODIFIED_TIMES_DIR/pending_recompute` (JSON) — written by save_modified_times
    when changes are detected; this script consumes it.

Flow:
  1. If marker is absent, exit 0 (nothing to do — most invocations).
  2. If marker is present, read it (logged for traceability), run the recompute
     pipeline via `ShroffCelegansModels.run_recompute_pipeline()`.
  3. On success, delete the marker. On failure, leave the marker so the next
     invocation retries.
"""

using Dates: now, format, Dates
using ShroffCelegansModels
using ShroffCelegansModels.JSON3

const _scripts_dir = joinpath(@__DIR__, "..", "..", "scripts")
include(joinpath(_scripts_dir, "export_meshscatter_static.jl"))
include(joinpath(_scripts_dir, "generate_pipeline_movie.jl"))
include(joinpath(_scripts_dir, "combined_meshscatter_movie.jl"))

const MARKER_NAME = "pending_recompute"

function read_marker(path::AbstractString)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

# Trigger a recompute when annotation_changes.h5 (the edits written by the web
# service's persist server) is newer than the last pipeline run. The marker
# mechanism only covers /nearline dataset (annotation CSV + lattice) mtime
# changes, so without this an annotation edit would never trigger a recompute.
# "Last run" = newest mtime among the top-level files in the recompute output
# dir (the pipeline rewrites those each run; archived/checkpoint subdirs skipped).
function annotation_changes_pending()
    changes = get(ENV, "ANNOTATION_CHANGES_PATH", "/data/annotations/annotation_changes.h5")
    isfile(changes) || return false
    output_dir = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute")
    isdir(output_dir) || return true   # never run before
    outputs = filter(isfile, joinpath.(output_dir, readdir(output_dir)))
    isempty(outputs) && return true
    return mtime(changes) > maximum(mtime, outputs)
end

function _generate_pipeline_visualizations(
    h5_path::AbstractString,
    unsmoothed_h5_path::Union{Nothing, AbstractString} = nothing,
)
    output_dir = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute")
    pipeline_run = Dates.format(Dates.unix2datetime(mtime(h5_path)), "yyyy-mm-ddTHH:MM:SS")
    avg_dict = ShroffCelegansModels.load_latest_average_annotations(
        default_filename = basename(h5_path),
        dir = dirname(h5_path),
    )

    export_path = joinpath(output_dir, "meshscatter_latest.html")
    export_meshscatter_static(avg_dict; output_path = export_path)

    for view in (:yz, :xz)
        movie_path = joinpath(output_dir, "movie_$(view).mp4")
        generate_meshscatter_movie(avg_dict; output_path = movie_path, view)
        @info "Movie written" view movie_path
    end

    # Combined pre-twitch + post-twitch movie (spans mpfc 20→751) in the same
    # meshscatter style. Pre-twitch frames are prepended and aligned to the
    # post-twitch H2 seam-cell midpoint; see combined_meshscatter_movie.jl.
    for view in (:yz, :xz)
        movie_path = joinpath(output_dir, "combined_movie_$(view).mp4")
        try
            generate_combined_meshscatter_movie(avg_dict; output_path = movie_path, view)
            @info "Combined movie written" view movie_path
        catch err
            @warn "Combined movie generation failed (other outputs still valid)" view err
        end
    end

    # Unsmoothed track: same movies, generated from the unsmoothed HDF5, for
    # QA/debug comparison against the smoothed movies above. Best-effort — a
    # failure here must never affect the smoothed outputs.
    if unsmoothed_h5_path !== nothing
        try
            unsmoothed_avg_dict = ShroffCelegansModels.load_latest_average_annotations(
                default_filename = basename(unsmoothed_h5_path),
                dir = dirname(unsmoothed_h5_path),
            )
            for view in (:yz, :xz)
                movie_path = joinpath(output_dir, "unsmoothed_movie_$(view).mp4")
                generate_meshscatter_movie(unsmoothed_avg_dict; output_path = movie_path, view)
                @info "Unsmoothed movie written" view movie_path
            end
            for view in (:yz, :xz)
                movie_path = joinpath(output_dir, "unsmoothed_combined_movie_$(view).mp4")
                generate_combined_meshscatter_movie(unsmoothed_avg_dict; output_path = movie_path, view)
                @info "Unsmoothed combined movie written" view movie_path
            end
        catch err
            @warn "Unsmoothed movie generation failed (smoothed outputs still valid)" err
        end
    end

    movie_created = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    status = (; pipeline_run, movie_created)
    open(joinpath(output_dir, "pipeline_status.json"), "w") do io
        JSON3.write(io, status)
    end

    # The pipeline's 8c step already wrote index.html, but that ran *before*
    # the movies/meshscatter HTML above were (re)generated, so the static
    # /recompute/ listing captured the previous run's mtimes for those files.
    # Regenerate it here, after the visualizations exist, so the page reflects
    # the freshly-written movie/HTML mtimes.
    try
        ShroffCelegansModels._write_recompute_index(output_dir)
        @info "Regenerated recompute index after visualizations" output_dir
    catch err
        @warn "Post-visualization index regeneration failed (outputs still valid)" err
    end

    @info "Pipeline visualizations complete" export_path pipeline_run movie_created
end

function run_pipeline(marker)
    kinds = if haskey(marker, :kinds)
        String[String(k) for k in marker[:kinds]]
    else
        # Legacy markers from the prior slice didn't carry `kinds`. Assume both.
        ["annotation", "lattice"]
    end
    # N_TIMEPOINTS env override lets the test CronJob run a faster validation
    # (e.g. N_TIMEPOINTS=51) without changing source defaults. Production uses
    # 371 (or whatever the deployment sets).
    n_timepoints = parse(Int, get(ENV, "N_TIMEPOINTS", "371"))
    result = ShroffCelegansModels.run_recompute_pipeline(; kinds = kinds, n_timepoints)
    @info "Pipeline produced artifacts" result

    try
        _generate_pipeline_visualizations(result.h5_path, result.unsmoothed_h5_path)
    catch err
        @warn "Visualization generation failed (pipeline outputs still valid)" err
    end

    return true
end

function main()
    output_dir = get(ENV, "MODIFIED_TIMES_DIR", "/data/annotations/modified_times")
    marker_path = joinpath(output_dir, MARKER_NAME)

    marker = read_marker(marker_path)
    if marker === nothing && annotation_changes_pending()
        @info "annotation_changes.h5 is newer than the last recompute output — triggering recompute"
        marker = Dict(:kinds => ["annotation", "lattice"], :reason => "annotation_changes_newer")
    end
    if marker === nothing
        @info "No recompute marker and annotation_changes.h5 not newer — nothing to do" marker_path
        return
    end

    @info "Running recompute pipeline" marker_path marker
    started_at = now()
    success = try
        run_pipeline(marker)
    catch err
        @error "Recompute pipeline raised — leaving marker in place for retry" err
        false
    end
    elapsed = now() - started_at

    if success
        @info "Pipeline succeeded — clearing marker" elapsed
        rm(marker_path; force=true)
    else
        @warn "Pipeline did not succeed — marker retained for next run" elapsed
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
