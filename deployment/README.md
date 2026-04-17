# shroff-data — Container Build & Deployment Package

This directory contains everything needed to rebuild the shroff-data container image
and deploy it on OpenShift. The `ShroffCelegansModels.jl/` source here is the **patched version**
already running in production — see "What was changed" below.

## Directory contents

```
shroff-data-for-mark/
├── README.md                  # This file
├── Dockerfile                 # Builds the Julia image
├── ShroffCelegansModels.jl/   # Patched source (ready to build from)
├── shroff-data/               # K8s manifests — production namespace (shroff-data.int.janelia.org)
│   ├── 00-namespace.yaml
│   ├── 01-configmap.yaml
│   ├── 02-pvc.yaml
│   ├── 03-deployment.yaml
│   ├── 04-route.yaml
│   ├── 05-buildconfig.yaml    ← BuildConfig is also here, inside the full manifest set
│   └── 06-scc.yaml
└── shroff-data-test/          # K8s manifests — test namespace (shroff-data-test.int.janelia.org)
    ├── 00-namespace.yaml
    ├── 01-configmap.yaml
    ├── 02-pvc.yaml
    ├── 03-deployment.yaml
    ├── 04-route.yaml
    ├── 05-buildconfig.yaml
    └── 06-scc.yaml
```

## Deploying from scratch

### Production namespace (`shroff-data`)

```bash
oc apply -f shroff-data/00-namespace.yaml
oc apply -f shroff-data/06-scc.yaml
oc apply -f shroff-data/05-buildconfig.yaml
oc start-build shroff-data --from-dir=. -n shroff-data --follow
oc apply -f shroff-data/01-configmap.yaml
oc apply -f shroff-data/02-pvc.yaml
oc apply -f shroff-data/03-deployment.yaml
oc apply -f shroff-data/04-route.yaml
```

### Test namespace (`shroff-data-test`)

```bash
oc apply -f shroff-data-test/00-namespace.yaml
oc apply -f shroff-data-test/06-scc.yaml
oc apply -f shroff-data-test/05-buildconfig.yaml
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow
oc apply -f shroff-data-test/01-configmap.yaml
oc apply -f shroff-data-test/02-pvc.yaml
oc apply -f shroff-data-test/03-deployment.yaml
oc apply -f shroff-data-test/04-route.yaml
```

The build takes ~15-30 minutes (Julia precompiles all packages).
The running pod restarts automatically once the new image is pushed.

## Rebuilding after a code update

Update the source in `ShroffCelegansModels.jl/`, then trigger a new build in whichever namespace you want to update:

```bash
# Production
oc start-build shroff-data --from-dir=. -n shroff-data --follow

# Test
oc start-build shroff-data-test --from-dir=. -n shroff-data-test --follow
```

## URLs

| Environment | URL |
|-------------|-----|
| Production  | https://shroff-data.int.janelia.org |
| Test        | https://shroff-data-test.int.janelia.org |

---

## What was changed from the original source

Four changes were needed to make the code work in a container.
**Please consider merging these back into your repo** so future builds are clean.

---

### 1. `src/demo_averaging/loading.jl` — config_path fallback

Your code only defines `config_path` for hostnames `"vm7249"` and `"KITTISOPIKULM-2"`.
A container has a random pod name, so `config_path` was never set and the app crashed.

Added an `else` branch using `@__DIR__` (works anywhere, not hostname-dependent):

```julia
# Add this else branch:
else
    const config_path = joinpath(@__DIR__, "..", "..", "config_2024_09_05_v1.json")
end
```

---

### 2. `scripts/launch_show_average_annotations.jl` — NFS alias fallback

`alias_cache_unix("/nearline/shroff")` was only called for hostname `"vm7249"`.
Added an `else` branch so it always runs in containers:

```julia
# Add this else branch:
else
    alias_cache_unix("/nearline/shroff")
end
```

---

### 3. `web/scripts/*.jl` — bind address (all 7 scripts)

Scripts were binding the HTTP server to `Sockets.getaddrinfo("shroff-data.int.janelia.org")`.
In a container that resolves to the load balancer IP (not a local interface), so the server crashed.

Change to `"0.0.0.0"` in all 7 scripts:

```julia
# Before:
Server(app, Sockets.getaddrinfo("shroff-data.int.janelia.org") |> string, 8180; ...)
# After:
Server(app, "0.0.0.0", 8180; ...)
```

---

### 4. `web/scripts/*.jl` — WebSocket proxy URL (all 7 scripts)

Bonito's `proxy_url` was hardcoded to `https://shroff-data.int.janelia.org/...`.
This tells the browser where to connect for WebSockets — it needs to match whatever
hostname the user is actually visiting (test route vs production).

Change to read from an environment variable:

```julia
# Before:
proxy_url="https://shroff-data.int.janelia.org/show_average_annotations/"
# After:
proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/show_average_annotations/"
```

The container sets `SHROFF_HOST` automatically from the deployment manifest. You don't need to set it yourself.

---

## Accessing HDF5 data files on the server

The HDF5 data files (previously at `/var/www/shroff/data/` on vm7249) are stored on a
persistent volume mounted at `/var/www/shroff/data` inside the nginx container.

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
