# Triggering the recompute pipeline manually

The recompute pipeline normally runs via the `recompute-on-mtime-change`
CronJob (daily at 03:00 ET, see `deployment/shroff-data-test/09-cronjob-recompute.yaml`).
The CronJob runs `web/scripts/run_recompute_if_needed.jl`, which is a
**conditional** entry point: it only executes the pipeline when one of these
is true, otherwise it exits 0 and does nothing —

1. a `pending_recompute` marker file exists at
   `/data/annotations/modified_times/pending_recompute`, **or**
2. `/data/annotations/annotation_changes.h5` is newer than the newest
   top-level file in `/data/annotations/recompute/` (i.e. an annotation edit
   landed since the last run).

So creating a Job from the CronJob is **not** enough on its own — if neither
condition holds, the job spins up Julia and immediately exits. To force a run
you must first plant the marker.

## Step-by-step (test namespace)

`oc` is provided by the `deployment/` pixi env (it is not on PATH). Examples
use the direct binary path; `cd deployment && pixi run oc …` works too.

```bash
OC=deployment/.pixi/envs/default/bin/oc
NS=shroff-data-test

# 1. Confirm login + project
$OC whoami
$OC project          # should be shroff-data-test

# 2. Write the pending_recompute marker.
#    The annotations PVC is mounted READ-ONLY in most containers of the web
#    deployment pod. Use a container that mounts it READ-WRITE — zscore-analysis
#    (or fix-ap-axis / debug-zscore) qualifies. The marker is JSON carrying the
#    `kinds` the pipeline should recompute.
POD=$($OC get pods -n $NS -l app=shroff-data-test -o name | head -1)
$OC exec -n $NS -c zscore-analysis $POD -- sh -c \
  'printf "%s" "{\"kinds\":[\"annotation\",\"lattice\"],\"reason\":\"manual_trigger\"}" \
   > /data/annotations/modified_times/pending_recompute'

# 3. Create a one-off Job from the CronJob (inherits N=371, --threads=16,
#    resources, volume mounts from the CronJob's jobTemplate).
$OC create job -n $NS --from=cronjob/recompute-on-mtime-change recompute-manual-$(date +%Y%m%d)

# 4. Watch it
$OC get pods -n $NS -l job-name=recompute-manual-$(date +%Y%m%d)
$OC logs -f -n $NS <pod-name>
```

## Notes

- **The script deletes the marker on success and leaves it on failure**, so a
  failed run is automatically retried by the next CronJob tick. If you abandon
  a manual run, remove the marker yourself or the next cron will re-trigger.
- A cold-cache full `N=371` run can take **hours** (the CronJob's
  `activeDeadlineSeconds` is 43200 = 12h). It archives the previous
  `/data/annotations/recompute/` outputs into a timestamped `archive_*` dir,
  then writes fresh averaged HDF5, CSVs, caches, the meshscatter HTML, and the
  yz/xz movies.
- For a faster validation run, override `N_TIMEPOINTS` (e.g. `51`). You can't
  set env on `oc create job --from`, so either edit the CronJob's env first or
  apply a one-off Job manifest with the override.
- `concurrencyPolicy: Forbid` only governs CronJob-spawned jobs; a manually
  created Job can still race a cron run on the shared
  `/data/annotations/recompute/checkpoint/` dir. Check for a running
  `recompute-on-mtime-change-*` pod before starting one (the cron schedule was
  lowered to once-daily specifically to avoid this — see task #41).
- **Production** (`shroff-data` namespace): substitute `NS=shroff-data` and the
  production PVC/CronJob names accordingly.
