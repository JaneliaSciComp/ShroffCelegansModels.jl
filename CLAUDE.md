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

### Rebuilding the container image (important)

These are **binary** BuildConfigs (`source.type: Binary`): the source is uploaded
from your working copy on every build. **Always pass `--from-dir=.` and run from the
repo root.** A build triggered without it (e.g. the web console "Start Build" button,
or a bare `oc start-build <name>`) receives no source archive and fails immediately
with `FetchSourceFailed` / "unable to extract binary build input".

```bash
# App image (~2 min) — routine code changes
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow

# Base image (~7 min) — only when Project.toml / Manifest.toml change
oc start-build shroff-data-test-base --from-dir=. -n shroff-data-test --follow
```

The repo is a Julia 1.12 workspace (`[workspace]` in the root `Project.toml`); all
members share the **root `Manifest.toml`** (now committed), which is what the app
build precompiles against. The deployment has no image trigger, so after a successful
build run `oc rollout restart deployment/shroff-data-test -n shroff-data-test` to pick
up the new `:latest` image.

### After deploying: run the endpoint check

After any rebuild/redeploy (`oc start-build …` + `oc rollout restart …`),
verify the web endpoints with:

```
julia --project=web deployment/post_deploy_check.jl [host]
```

It GETs every public Julia endpoint and checks both the HTTP status **and**
the response body for a Julia-exception signature — Bonito renders a
server-side error into the page with HTTP 200, so a status-only check misses
render-time failures. Host defaults to the test env; pass
`shroff-data.int.janelia.org` (or set `SHROFF_CHECK_HOST`) to check production.
Exit code is non-zero if any endpoint fails.

### Live data location (important)

The pipeline and web service read/write data on an **OpenShift PVC**, not the
local repo. The local `*.h5` files in the repo root (e.g.
`annotation_changes_*.h5`, `annotations_cache.h5`) are **stale downloaded
snapshots** — do not treat them as the source of truth. Check the PVC via `oc`
(e.g. `oc rsh` into the relevant pod or `oc exec`) for current state.

Key PVC paths (see `recompute_pipeline.jl` defaults):
- `/data/annotations/annotation_changes.h5` — edits written by the web service
- `/data/annotations/recompute/` — pipeline output dir (averaged HDF5, CSVs, caches)

## LSF (cluster job scheduler)

Janelia uses IBM LSF for batch jobs. Common queues: `gpu_l4` (NVIDIA L4 GPU nodes),
`local` (shared CPU nodes). Always bill to `-P scicompsoft`.

```bash
# Submit a job script
bsub < scripts/submit_movie_glmakie.bsub

# Check job status
bjobs <JOBID>

# Detailed job info / exit reason
bjobs -l <JOBID>

# Job history (after completion)
bhist -l <JOBID>
```

### GPU jobs (gpu_l4 queue)

Request a GPU with `-gpu "num=1:j_exclusive=yes"`. GLMakie requires a display
even for offscreen rendering — use `xvfb-run -a` to provide a virtual X11 display:

```bash
#BSUB -q gpu_l4
#BSUB -P scicompsoft
#BSUB -gpu "num=1:j_exclusive=yes"

xvfb-run -a julia --project=glmakie glmakie/scripts/generate_movie_glmakie.jl "$H5" "$OUT"
```

`JULIA_GLMAKIE_BACKEND=egl` alone is **not sufficient** — GLFW's `Init()` on Linux
defaults to `PLATFORM_X11` and fails without a display regardless of the EGL setting.

`JULIA_GLFW_PLATFORM=null` (a proposed fix in GLFW.jl issue #253) routes GLFW through
`PLATFORM_NULL` which uses OSMesa for off-screen rendering — but OSMesa is not installed
on these cluster nodes, so it fails with `OSMesa: Library not found`. The `xvfb-run -a`
workaround remains the only working option on Janelia HPC until OSMesa is available.

### Julia path on compute nodes

`julia` (from juliaup) may not be in PATH on compute nodes. Use the full path:

```bash
JULIA=/groups/scicompsoft/home/kittisopikulm/.juliaup/bin/julia
xvfb-run -a "$JULIA" --project=...
```

### Important: `/tmp` is node-local

Compute nodes do not share `/tmp` with the login node. Input files (h5, etc.) and
log output must be on a shared NFS path (e.g. `/groups/…` or `movies/`) — not `/tmp`.

```bash
#BSUB -o /groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl/movies/job_%J.log
#BSUB -e /groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl/movies/job_%J.err
```

### Reading job logs

Log files are written to the path specified by `-o`/`-e`. If those point to a shared
filesystem, read them directly. If they point to `/tmp` on the compute node, you
cannot SSH to compute nodes — use `bhist -l <JOBID>` for a summary, or always set
`-o`/`-e` to a shared path from the start.
