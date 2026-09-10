# Regenerating display/export artifacts only (no re-averaging)

When a change touches only the **display / export** code paths — e.g. the LR/DV
axis flip in `load_average_annotations` (`src/demo_averaging/average_annotations.jl`)
or the CSV sign/columns in `src/demo_averaging/explicit_export.jl` — the averaged
HDF5 (`edited_smoothed_average_annotations_*.h5`) is **numerically unchanged**.
Only the derived artifacts need rebuilding:

- `meshscatter_latest.html`
- `movie_yz.mp4`, `movie_xz.mp4`
- `combined_movie_yz.mp4`, `combined_movie_xz.mp4`
- `pretwitch_<ts>.csv`, `posttwitch_<ts>.csv`, `combined_<ts>.csv`
- `index.html`

A full recompute (`recompute-on-mtime-change`) would re-average every timepoint
(~hours). Instead run the **viz-only** regeneration (~minutes), which reuses the
newest averaged HDF5 in place: `web/scripts/regenerate_viz_only.jl`. It shares
the pipeline's own viz routine (`_generate_pipeline_visualizations` in
`web/scripts/run_recompute_if_needed.jl`) and steps 8/8a for the CSVs. The CSVs
are written with the timestamp parsed from the HDF5 name, so they **overwrite**
the existing trio rather than accumulating a new dated set. The averaged HDF5 and
the `avg_models_n<N>.h5` cache are read, never rewritten.

## Trigger (test namespace)

`oc` is provided by the `deployment/` pixi env (not on PATH).

```bash
OC=deployment/.pixi/envs/default/bin/oc
NS=shroff-data-test

# 1. Confirm login + project
$OC whoami
$OC project          # should be shroff-data-test

# 2. Make sure no recompute job is running (both write to
#    /data/annotations/recompute/).
$OC get pods -n $NS -l component=recompute-on-mtime-change

# 3. Apply the suspended CronJob template (idempotent), then create a Job from it.
$OC apply -f deployment/shroff-data-test/11-cronjob-regenerate-viz.yaml -n $NS
$OC create job -n $NS --from=cronjob/regenerate-viz regenerate-viz-$(date +%Y%m%d%H%M)

# 4. Watch it
$OC get pods -n $NS -l component=regenerate-viz -w
$OC logs -f -n $NS <pod-name>
```

## Notes

- The CronJob is **suspended** (`suspend: true`) — it never fires on a schedule;
  it's only a template for `oc create job --from`. No `pending_recompute` marker
  is involved (this bypasses the conditional recompute entry point and always
  regenerates).
- The Job mounts the annotations PVC **read-write** (most web-deployment
  containers mount it read-only). It runs the same `:latest` app image as the web
  deployment, so build + roll out the image *before* triggering.
- Toggles (env on the Job, or edit the CronJob): `REGEN_CSVS=0` to skip the CSV
  rebuild, `REGEN_VIZ=0` to skip the movies/HTML/index. `RECOMPUTE_OUTPUT_DIR`
  and `N_TIMEPOINTS` default to the production values.
- **Production** (`shroff-data`): a parallel manifest under `deployment/shroff-data/`
  can be added when needed; substitute `NS=shroff-data`.
