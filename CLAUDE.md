# ShroffCelegansModels.jl — agent notes

## OpenShift / deployment tooling

The OpenShift CLI (`oc`) is **not** on the system `PATH`. It is provided by a
**pixi environment in the `deployment/` folder** (`deployment/pixi.toml`
declares `openshift-cli`).

To run `oc` against the cluster, either:

- `cd deployment && pixi run oc <args>`, or
- call the binary directly: `deployment/.pixi/envs/default/bin/oc <args>`

The production and test deployments live under:

- `deployment/shroff-data/` — production
- `deployment/shroff-data-test/` — test environment (`shroff-data-test.int.janelia.org`)

### Live data location (important)

The pipeline and web service read/write data on an **OpenShift PVC**, not the
local repo. The local `*.h5` files in the repo root (e.g.
`annotation_changes_*.h5`, `annotations_cache.h5`) are **stale downloaded
snapshots** — do not treat them as the source of truth. Check the PVC via `oc`
(e.g. `oc rsh` into the relevant pod or `oc exec`) for current state.

Key PVC paths (see `recompute_pipeline.jl` defaults):
- `/data/annotations/annotation_changes.h5` — edits written by the web service
- `/data/annotations/recompute/` — pipeline output dir (averaged HDF5, CSVs, caches)
