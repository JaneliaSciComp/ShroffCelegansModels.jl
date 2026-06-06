# shroff-data — Container Build & Deployment Package

This directory contains everything needed to rebuild the shroff-data container image
and deploy it on OpenShift. The source is the `ShroffCelegansModels.jl` repo itself —
all patches previously applied in the separate build package have been merged back in.

## Directory contents

```
ShroffCelegansModels.jl/           # Repo root — also the Docker build context
├── .dockerignore                  # Excludes deployment/, backup_*/, cache_*/, services/
├── deployment/
│   ├── README.md                  # This file
│   ├── Dockerfile                 # App image — builds FROM the base image (~2 min)
│   ├── Dockerfile.base            # Base image — installs deps & precompiles external packages (~7 min)
│   ├── shroff-data/               # K8s manifests — production namespace (shroff-data.int.janelia.org)
│   │   ├── 00-namespace.yaml
│   │   ├── 01-configmap.yaml
│   │   ├── 02-pvc.yaml
│   │   ├── 03-deployment.yaml
│   │   ├── 04-route.yaml
│   │   ├── 05-buildconfig.yaml
│   │   └── 06-scc.yaml
│   └── shroff-data-test/          # K8s manifests — test namespace (shroff-data-test.int.janelia.org)
│       ├── 00-namespace.yaml
│       ├── 01-configmap.yaml
│       ├── 02-pvc.yaml
│       ├── 03-deployment.yaml
│       ├── 04-route.yaml
│       ├── 05-buildconfig.yaml
│       ├── 05b-buildconfig-base.yaml
│       └── 06-scc.yaml
```

All `oc` commands below must be run from the **repo root** (`ShroffCelegansModels.jl/`),
which is the Docker build context.

## Prerequisites: `oc` via pixi

The OpenShift CLI is provided by [pixi](https://pixi.sh/) — see `deployment/pixi.toml`
(`openshift-cli` from conda-forge). The system `PATH` does not include `oc`; invoke it
through pixi instead.

Pick one of the following:

```bash
# Option A — prefix every command (works from repo root)
pixi run --manifest-path deployment/pixi.toml oc whoami

# Option B — drop into a shell with oc on PATH
cd deployment && pixi shell    # then `oc ...` works; `cd ..` back to repo root for builds

# Option C — alias in your own shell
alias oc='pixi run --manifest-path '"$PWD"'/deployment/pixi.toml oc'
```

The `oc` invocations in the rest of this README assume one of the above is in effect.
Log in with `oc login ...` against the cluster before continuing.

## Two-stage build

The build is split into two images to keep routine deploys fast:

| Image | Dockerfile | Build time | Rebuild when |
|-------|-----------|------------|--------------|
| `shroff-data-test-base` | `Dockerfile.base` | ~7 min | `Project.toml` or `Manifest.toml` changes |
| `shroff-data-test` | `Dockerfile` | ~2 min | Any code change |

`Dockerfile.base` installs system packages and precompiles all external Julia packages
using a stub for the local `ShroffCelegansModels` package. `Dockerfile` builds `FROM`
that base, copies the real source, and re-precompiles only the local package.

## Deploying from scratch

### Test namespace (`shroff-data-test`)

```bash
oc apply -f deployment/shroff-data-test/00-namespace.yaml
oc apply -f deployment/shroff-data-test/06-scc.yaml

# Build the base image first (~7 min)
oc apply -f deployment/shroff-data-test/05b-buildconfig-base.yaml
oc start-build shroff-data-test-base --from-dir=. -n shroff-data-test --follow

# Then build the app image (~2 min)
oc apply -f deployment/shroff-data-test/05-buildconfig.yaml
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow

oc apply -f deployment/shroff-data-test/01-configmap.yaml
oc apply -f deployment/shroff-data-test/02-pvc.yaml
oc apply -f deployment/shroff-data-test/03-deployment.yaml
oc apply -f deployment/shroff-data-test/04-route.yaml
```

### Production namespace (`shroff-data`)

> **Note:** The production namespace does not yet have a base image BuildConfig.
> Use the single-stage build until it is set up (see test namespace as a reference).

```bash
oc apply -f deployment/shroff-data/00-namespace.yaml
oc apply -f deployment/shroff-data/06-scc.yaml
oc apply -f deployment/shroff-data/05-buildconfig.yaml
oc start-build shroff-data --from-dir=. -n shroff-data --follow
oc apply -f deployment/shroff-data/01-configmap.yaml
oc apply -f deployment/shroff-data/02-pvc.yaml
oc apply -f deployment/shroff-data/03-deployment.yaml
oc apply -f deployment/shroff-data/04-route.yaml
```

## Rebuilding after a code update (~2 min)

For routine code changes, only the app image needs to be rebuilt:

```bash
# Test
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow

# Production
oc start-build shroff-data --from-dir=. -n shroff-data --follow
```

## Rebuilding after a dependency update (~7 + 2 min)

When `Project.toml` or `Manifest.toml` changes, rebuild the base first, then the app:

```bash
# Test
oc start-build shroff-data-test-base --from-dir=. -n shroff-data-test --follow
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow
```

## URLs

| Environment | URL |
|-------------|-----|
| Production  | https://shroff-data.int.janelia.org |
| Test        | https://shroff-data-test.int.janelia.org |

---

## Patches merged into the repo

These four changes were previously applied to a separate build copy and have since been
merged back into the repo. They are noted here for reference.

### 1. `src/demo_averaging/loading.jl` — config_path fallback

Added an `else` branch so the config path resolves correctly in containers (which have
random pod names, not the expected hostnames):

```julia
else
    const config_path = joinpath(@__DIR__, "..", "..", "config", "linux", "config_2026_03_19_v2.json")
end
```

### 2. `scripts/launch_show_average_annotations.jl` — NFS alias fallback

Added an `else` branch so `alias_cache_unix("/nearline/shroff")` always runs in containers:

```julia
else
    alias_cache_unix("/nearline/shroff")
end
```

### 3. `web/scripts/*.jl` — bind address

Changed HTTP server bind address from `Sockets.getaddrinfo("shroff-data.int.janelia.org")`
to `"0.0.0.0"` so it binds to a local interface inside the container.

### 4. `web/scripts/*.jl` — WebSocket proxy URL

Changed Bonito's `proxy_url` from a hardcoded hostname to read from an environment variable,
allowing the same image to serve both production and test:

```julia
proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/show_average_annotations/"
```

`SHROFF_HOST` is set automatically by the deployment manifest.

---

## Accessing HDF5 data files on the server

The HDF5 data files are stored on a persistent volume mounted at `/var/www/shroff/data`
inside the nginx container.

### Get a shell to browse or edit files

```bash
# Production
oc exec -it -n shroff-data deployment/shroff-data -c nginx -- sh

# Test
oc exec -it -n shroff-data-test deployment/shroff-data-test -c nginx -- sh

# Files are at:
ls /var/www/shroff/data/
```

### Copy a file from your machine to the server

```bash
# Get the pod name first
oc get pods -n shroff-data

oc cp /local/path/to/file.h5 shroff-data/<pod-name>:/var/www/shroff/data/file.h5 -c nginx
```

### Copy a file from the server to your machine

```bash
oc cp shroff-data/<pod-name>:/var/www/shroff/data/file.h5 /local/path/to/file.h5 -c nginx
```

### Delete or rename files

```bash
oc exec -n shroff-data deployment/shroff-data -c nginx -- rm /var/www/shroff/data/old.h5
oc exec -n shroff-data deployment/shroff-data -c nginx -- mv /var/www/shroff/data/old.h5 /var/www/shroff/data/new.h5
```

### Notes

- Changes to the PVC are **persistent** — they survive pod restarts and redeployments
- The volume is shared between the nginx container (serves the files) and nothing else — Julia apps use `/nearline/shroff` for microscopy data, not this PVC
- `oc` CLI must be installed and you must be logged in to the OpenShift cluster
