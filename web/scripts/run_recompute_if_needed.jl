"""
Runs the recompute pipeline if the daily `save_modified_times.jl` CronJob
detected mtime changes since its previous snapshot.

Pattern: filesystem marker on the shared annotations PVC.

  - `MODIFIED_TIMES_DIR/pending_recompute` (JSON) — written by save_modified_times
    when changes are detected; this script consumes it.

Flow:
  1. If marker is absent, exit 0 (nothing to do — most invocations).
  2. If marker is present, read it (logged for traceability), run the recompute
     pipeline.
  3. On success, delete the marker. On failure, leave the marker so the next
     invocation retries.

The pipeline itself is currently a STUB (logs the marker contents and exits
successfully). Steps 1–6 of the recompute (get_avg_models → average_annotations
→ load_annotation_changes_cache → update_annotations_cache → write averaged
HDF5 → export DataFrame) land in a follow-up.
"""

using Dates: now, format
using ShroffCelegansModels.JSON3

const MARKER_NAME = "pending_recompute"

function read_marker(path::AbstractString)
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

# STUB. Replace this with the real pipeline (steps 1–6) when it's ready.
function run_pipeline(marker)
    @info "[STUB] Would run recompute pipeline" marker
    @info "[STUB] Pretending pipeline succeeded"
    return true
end

function main()
    output_dir = get(ENV, "MODIFIED_TIMES_DIR", "/data/annotations/modified_times")
    marker_path = joinpath(output_dir, MARKER_NAME)

    marker = read_marker(marker_path)
    if marker === nothing
        @info "No recompute marker — nothing to do" marker_path
        return
    end

    @info "Recompute marker found — running pipeline" marker_path marker
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
