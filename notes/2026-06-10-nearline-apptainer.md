# NEARLINE_BASE env var + Apptainer testing (2026-06-10)

## Problem

The LSF cluster nodes do not have `/nearline` mounted.  The pipeline reads C.
elegans CSV and lattice data from `/nearline/shroff/shrofflab/…`.  The 72 config
datasets have corresponding tar.gz archives already on the PVC at
`/data/annotations/csv_archives/`.

## Option A — NEARLINE_BASE env var (lightweight, no root needed)

Set `NEARLINE_BASE` to any local directory that mirrors the `/nearline` tree
before starting Julia.  All dataset paths and cache keys starting with
`/nearline` are transparently remapped at pipeline startup.

Two code sites read the env var:

| File | Location | Effect |
|------|----------|--------|
| `src/recompute_pipeline.jl` | `_NEARLINE_BASE` const + `_remap_nearline()` | normalises Windows cache keys **and** live Linux paths for staleness checks |
| `src/demo_averaging/read_config_json.jl` | `remap` lambda inside `read_config_json` | remaps `folder_path` before constructing each `NormalizedDataset` |

### Preparing the local mirror

```bash
# 1. Download the 72 archives from the PVC (run on a machine with oc access)
mkdir -p ~/nearline_local/shroff
for f in $(deployment/.pixi/envs/default/bin/oc exec -n shroff-data-test \
               deploy/shroff-data-test -- ls /data/annotations/csv_archives/); do
    deployment/.pixi/envs/default/bin/oc cp \
        shroff-data-test/$(oc get pod -l app=shroff-data-test -o name | head -1 | cut -d/ -f2):/data/annotations/csv_archives/$f \
        ~/nearline_local/shroff/$f
done

# 2. Unpack all archives into ~/nearline_local/shroff/
#    Archives were created with -C /nearline/shroff/ so they contain
#    relative paths like  shrofflab/OD1599_NU/…
cd ~/nearline_local/shroff
for f in *.tar.gz; do tar xzf "$f"; done
# Result: ~/nearline_local/shroff/shrofflab/OD1599_NU/…/RegB/…
```

### Running the pipeline with NEARLINE_BASE

```bash
NEARLINE_BASE=~/nearline_local/shroff \
    julia --project=web web/scripts/run_recompute_if_needed.jl
```

Or when calling `run_recompute_pipeline` directly from Julia:

```julia
ENV["NEARLINE_BASE"] = expanduser("~/nearline_local/shroff")
# must be set BEFORE `using ShroffCelegansModels` (const is read at load time)
using ShroffCelegansModels
ShroffCelegansModels.run_recompute_pipeline()
```

## Option B — Apptainer container (simulates /nearline, zero code changes needed)

Useful for testing code that predates the NEARLINE_BASE env var, or for
validating that NEARLINE_BASE works correctly.  Apptainer 1.5.0 is available at
`/usr/bin/apptainer`.

```bash
# Bind the local mirror to the canonical /nearline/shroff path
apptainer exec \
    --bind ~/nearline_local/shroff:/nearline/shroff \
    --bind /groups:/groups \
    /path/to/julia.sif \
    julia --project=/groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl/web \
          /groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl/web/scripts/run_recompute_if_needed.jl
```

If you don't have a Julia SIF image, build one:

```bash
# Minimal SIF wrapping the official Julia tarball
cat > julia.def <<'EOF'
Bootstrap: docker
From: julia:1.12
EOF
apptainer build julia.sif julia.def
```

With this bind mount, the container sees `/nearline/shroff/shrofflab/…` exactly as the
production OpenShift pod does — no code changes required, good for integration testing.

## LSF parallelisation (future)

With NEARLINE_BASE in place, separate LSF jobs can process subsets of the 72
datasets concurrently (each job sets `NEARLINE_BASE` and a config subset) and
write their output to separate directories on `/groups/scicompsoft/`.  See git
log for commit `067efb8` ("Parallelize across datasets") for a prior sketch.
