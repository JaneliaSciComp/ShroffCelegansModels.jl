#!/usr/bin/env bash
#
# Monitor a recompute pipeline run on shroff-data-test.
#
# Usage:
#   deployment/monitor_recompute.sh [JOB_NAME] [REFRESH_SECONDS]
#
#   JOB_NAME         defaults to the most recently created recompute-manual-* job
#   REFRESH_SECONDS  poll interval (default 30)
#
# The recompute pipeline logs progress as [N/8] step markers and "Step done"
# timing lines; this script surfaces those plus job/pod status, and exits when
# the job Completes or Fails.

set -euo pipefail

NS="${NAMESPACE:-shroff-data-test}"
JOB="${1:-}"
REFRESH="${2:-30}"

# Resolve `oc` — prefer the pixi env binary (no per-call pixi overhead),
# fall back to `pixi run oc` from the deployment dir.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OC_BIN="$SCRIPT_DIR/.pixi/envs/default/bin/oc"
if [[ -x "$OC_BIN" ]]; then
  oc() { "$OC_BIN" "$@"; }
else
  oc() { ( cd "$SCRIPT_DIR" && pixi run oc "$@" ); }
fi

# Default to the newest recompute-manual-* job.
if [[ -z "$JOB" ]]; then
  JOB="$(oc get jobs -n "$NS" \
          --sort-by=.metadata.creationTimestamp \
          -o jsonpath='{range .items[*]}{.metadata.name}{"\n"}{end}' \
        | grep '^recompute-manual-' | tail -1 || true)"
  if [[ -z "$JOB" ]]; then
    echo "No recompute-manual-* job found in namespace $NS." >&2
    exit 1
  fi
  echo "No job specified — monitoring latest: $JOB"
fi

while true; do
  POD="$(oc get pods -n "$NS" -l "job-name=$JOB" \
          -o jsonpath='{.items[0].metadata.name}' 2>/dev/null || true)"

  clear
  echo "=========================================================="
  echo " Recompute monitor — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " namespace=$NS  job=$JOB"
  echo "=========================================================="

  echo
  echo "── Job ─────────────────────────────────────────────────"
  oc get job "$JOB" -n "$NS" \
    -o custom-columns='STATUS:.status.conditions[0].type,COMPLETIONS:.status.succeeded,FAILED:.status.failed,START:.status.startTime' \
    2>/dev/null || echo "  (job not found yet)"

  echo
  echo "── Pod ─────────────────────────────────────────────────"
  if [[ -n "$POD" ]]; then
    oc get pod "$POD" -n "$NS" \
      -o custom-columns='NAME:.metadata.name,READY:.status.containerStatuses[0].ready,PHASE:.status.phase,RESTARTS:.status.containerStatuses[0].restartCount' \
      2>/dev/null || echo "  (pod gone)"
  else
    echo "  (no pod yet)"
  fi

  if [[ -n "$POD" ]]; then
    echo
    echo "── Pipeline step (latest [N/8]) ────────────────────────"
    oc logs "$POD" -n "$NS" 2>/dev/null \
      | grep -E '\[[0-9]+(\.[0-9]+)?[a-z]?/8\]' | tail -3 || echo "  (no step markers yet)"

    echo
    echo "── Completed steps (Step done) ─────────────────────────"
    oc logs "$POD" -n "$NS" 2>/dev/null \
      | grep -E 'Step done' | tail -6 || echo "  (none yet)"

    echo
    echo "── Last log lines ──────────────────────────────────────"
    oc logs "$POD" -n "$NS" --tail=8 2>/dev/null || echo "  (no logs yet)"
  fi

  # Exit on terminal job state.
  STATE="$(oc get job "$JOB" -n "$NS" \
            -o jsonpath='{.status.conditions[0].type}' 2>/dev/null || true)"
  if [[ "$STATE" == "Complete" ]]; then
    echo
    echo "✅ Job $JOB COMPLETED."
    exit 0
  elif [[ "$STATE" == "Failed" ]]; then
    echo
    echo "❌ Job $JOB FAILED."
    exit 1
  fi

  echo
  echo "(refreshing in ${REFRESH}s — Ctrl-C to stop; the job keeps running)"
  sleep "$REFRESH"
done
