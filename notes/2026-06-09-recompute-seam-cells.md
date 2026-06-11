# 2026-06-09 — Recompute run (prod annotation edits + seam cells)

## What's running
Manual recompute Job on **shroff-data-test**, triggered after:
- Uploading prod `annotation_changes.h5` (last edit: `RIGR` @ tp 007,
  `CND-1_RedUntwisting_A/.../Pos3/.../RegB`, 2026-06-08 14:47 UTC) over the test PVC.
- Adding the seam-cell injection (`avg_dict["seam_cells"] = seam_cells_as_annotations(avg_models)`)
  before smoothing — commit `90a5206`, image build #55.

Job name: `recompute-manual-20260609-070036`
(created via `oc create job --from=cronjob/recompute-on-mtime-change`).
Trigger confirmed in logs: `annotation_changes.h5 is newer ... reason => annotation_changes_newer`.

## Check status
`oc` is via pixi in `deployment/`. Run from the deployment dir (or use `cd deployment && pixi run ...`):

```bash
cd deployment

# Pod for this job (name may change if it restarted):
pixi run oc get pods -n shroff-data-test -l job-name=recompute-manual-20260609-070036

# Tail the pipeline log:
pixi run oc logs recompute-manual-20260609-070036-5tngt -n shroff-data-test | tail

# Job done?
pixi run oc get job recompute-manual-20260609-070036 -n shroff-data-test
```

## Confirm fresh output (with seam cells)
List newest averaged-annotation HDF5 on the PVC (rsh into any pod that mounts it):

```bash
pixi run oc rsh deploy/<pod-with-annotations-pvc> \
  ls -lt /data/annotations/recompute/edited_smoothed_average_annotations_*.h5
```

The newest file should be dated 2026-06-09 and, when loaded, contain a
`seam_cells` group alongside the per-experiment annotation groups. The
meshscatter web app (`load_latest_average_annotations`) picks it up on the
next service restart.
