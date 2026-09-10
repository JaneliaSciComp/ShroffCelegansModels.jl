"""
Regenerate the recompute pipeline's **display/export artifacts only**, from the
already-averaged HDF5 — **without** re-running the expensive averaging.

Use this whenever a fix changes only the display/export code paths (e.g. the
LR/DV axis flip in `load_average_annotations` / `explicit_export.jl`) and the
averaged HDF5 itself is unchanged. A full `run_recompute_if_needed.jl` run would
re-average all timepoints (~hours); this reuses the newest
`edited_smoothed_average_annotations_*.h5` in place and rebuilds:

  - `meshscatter_latest.html`                      (export_meshscatter_static)
  - `movie_yz.mp4`, `movie_xz.mp4`                 (generate_meshscatter_movie)
  - `combined_movie_yz.mp4`, `combined_movie_xz.mp4` (generate_combined_meshscatter_movie)
  - `pretwitch_<ts>.csv`, `posttwitch_<ts>.csv`, `combined_<ts>.csv`
      (resave_for_ben → write_combined_explicit_csvs)
  - `index.html`                                   (_write_recompute_index)

When the run also produced the unsmoothed track (newest
`edited_unsmoothed_average_annotations_*.h5`), its counterparts are rebuilt too:

  - `unsmoothed_movie_{yz,xz}.mp4`, `unsmoothed_combined_movie_{yz,xz}.mp4`
  - `unsmoothed_{pretwitch,posttwitch,combined}_<ts>.csv`

Set `REGEN_VIZ=0 REGEN_CSVS=0` to rebuild **only** the unsmoothed movies —
useful after a fix that affected just that track (each 371-frame movie takes
~20 min, so skipping the smoothed set roughly halves the run).

The CSVs are written with the timestamp parsed from the HDF5 filename, so they
**overwrite** the existing trio (same names) rather than piling up a new dated
set. The averaged HDF5 and the `avg_models_n<N>.h5` cache are read, never
rewritten. The viz half is shared verbatim with the pipeline by reusing
`_generate_pipeline_visualizations` from `run_recompute_if_needed.jl`.

Usage:
    julia --project=web web/scripts/regenerate_viz_only.jl

Environment:
    RECOMPUTE_OUTPUT_DIR   output/HDF5 dir      (default /data/annotations/recompute)
    N_TIMEPOINTS           avg_models N to load (default 371)  — only used to log;
                           the newest avg_models_n*.h5 is loaded by mtime.
    REGEN_CSVS             "0" to skip the CSV rebuild   (default on)
    REGEN_VIZ             "0" to skip movies/HTML/index  (default on)
    REGEN_UNSMOOTHED      "0" to skip the unsmoothed track entirely (default on)
"""

using ShroffCelegansModels

# Reuse the pipeline's viz routine (movies + meshscatter HTML + index) and its
# script-level includes. run_recompute_if_needed.jl guards its own `main()`
# behind `abspath(PROGRAM_FILE) == @__FILE__`, so including it here defines
# `_generate_pipeline_visualizations` (and pulls in the movie/HTML scripts)
# without running the marker-driven pipeline.
include(joinpath(@__DIR__, "run_recompute_if_needed.jl"))

# Matches the `yyyy_mm_dd_HHMMSS` stamp the pipeline embeds in every dated output.
const _TS_RE = r"(\d{4}_\d{2}_\d{2}_\d{6})"

# The two averaged-annotation tracks the pipeline writes per run. Neither
# prefix matches the other's filenames, so each glob sees only its own track.
const _SMOOTHED_PREFIX = "edited_smoothed_average_annotations_"
const _UNSMOOTHED_PREFIX = "edited_unsmoothed_average_annotations_"

"Newest `<prefix>*.h5` in `output_dir` by mtime, or `nothing` if none match."
function _latest_h5(output_dir::AbstractString, prefix::AbstractString)
    isdir(output_dir) || return nothing
    candidates = String[
        joinpath(output_dir, f) for f in readdir(output_dir)
        if startswith(f, prefix) && endswith(f, ".h5") && !occursin(".tmp.", f)
    ]
    return isempty(candidates) ? nothing : argmax(mtime, candidates)
end

"Newest smoothed averaged HDF5 in `output_dir`; errors if there is none."
function _latest_average_h5(output_dir::AbstractString)
    isdir(output_dir) || error("regenerate_viz_only: output dir not found: $output_dir")
    path = _latest_h5(output_dir, _SMOOTHED_PREFIX)
    path === nothing &&
        error("regenerate_viz_only: no $(_SMOOTHED_PREFIX)*.h5 in $output_dir")
    return path
end

"""
    _regenerate_csvs(h5_path, output_dir; prefix = "")

Rebuild the explicit `pretwitch/posttwitch/combined` CSV trio from `h5_path`,
mirroring steps 8/8a of `run_recompute_pipeline` (resave_for_ben →
write_combined_explicit_csvs, then remove the intermediates). The `date_str` is
taken from the HDF5 filename so the trio overwrites the existing files.

Pass `prefix = "unsmoothed"` alongside the unsmoothed HDF5 to rebuild that
track's trio (`unsmoothed_pretwitch_<ts>.csv`, …).
"""
function _regenerate_csvs(
    h5_path::AbstractString,
    output_dir::AbstractString;
    prefix::AbstractString = "",
)
    m = match(_TS_RE, basename(h5_path))
    ts = m === nothing ? Dates.format(Dates.now(), "yyyy_mm_dd_HHMMSS") : m.captures[1]

    ben_csv = replace(h5_path, ".h5" => "_for_ben.csv")
    @info "[csv] Building intermediate _for_ben CSV" ben_csv
    ShroffCelegansModels.resave_for_ben(h5_path;
        target_filename = ben_csv, time_range = (381, 751),
        canonicalize_cell_names = true)

    avg_models = ShroffCelegansModels.load_latest_avg_models(dir = output_dir)
    @info "[csv] Writing explicit CSVs (overwriting the ts-matched trio)" ts prefix n_models=length(avg_models)
    res = ShroffCelegansModels.write_combined_explicit_csvs(;
        output_dir = output_dir,
        avg_models = avg_models,
        ben_csv_path = ben_csv,
        prefix = prefix,
        date_str = ts,
        add_unsmoothed_seam_cells = false,
    )
    @info "[csv] Wrote explicit CSVs" res.pretwitch_path res.posttwitch_path res.combined_path

    for f in (ben_csv,
              replace(ben_csv, ".csv" => "_ryan_duplicates.csv"),
              replace(ben_csv, ".csv" => "_ryan_stats.csv"))
        isfile(f) && rm(f; force = true)
    end
    return res
end

function main()
    output_dir = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute")
    n_timepoints = parse(Int, get(ENV, "N_TIMEPOINTS", "371"))
    do_csvs = get(ENV, "REGEN_CSVS", "1") != "0"
    do_viz  = get(ENV, "REGEN_VIZ", "1") != "0"
    do_unsmoothed = get(ENV, "REGEN_UNSMOOTHED", "1") != "0"

    h5_path = _latest_average_h5(output_dir)
    unsmoothed_h5_path = do_unsmoothed ? _latest_h5(output_dir, _UNSMOOTHED_PREFIX) : nothing
    if do_unsmoothed && unsmoothed_h5_path === nothing
        @warn "No unsmoothed HDF5 found — skipping the unsmoothed track" output_dir prefix=_UNSMOOTHED_PREFIX
    end
    @info "regenerate_viz_only starting" output_dir h5_path unsmoothed_h5_path n_timepoints do_csvs do_viz do_unsmoothed

    if do_csvs
        _regenerate_csvs(h5_path, output_dir)
        unsmoothed_h5_path === nothing ||
            _regenerate_csvs(unsmoothed_h5_path, output_dir; prefix = "unsmoothed")
    end

    if do_viz
        # Regenerates movies (both tracks) + meshscatter HTML, then rewrites
        # index.html last so it reflects the fresh CSV *and* movie mtimes.
        _generate_pipeline_visualizations(h5_path, unsmoothed_h5_path)
    else
        # Smoothed viz pass skipped: still rebuild the unsmoothed movies if asked
        # (REGEN_VIZ=0 REGEN_CSVS=0 regenerates *only* that set), then refresh
        # the index for whatever did change.
        if unsmoothed_h5_path !== nothing
            try
                _generate_unsmoothed_movies(unsmoothed_h5_path, output_dir)
            catch err
                @warn "Unsmoothed movie generation failed (other outputs still valid)" err
            end
        end
        if do_csvs || unsmoothed_h5_path !== nothing
            try
                ShroffCelegansModels._write_recompute_index(output_dir)
                @info "Regenerated recompute index (no smoothed viz pass)" output_dir
            catch err
                @warn "Index regeneration failed (other outputs still valid)" err
            end
        end
    end

    @info "regenerate_viz_only complete" output_dir h5_path unsmoothed_h5_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
