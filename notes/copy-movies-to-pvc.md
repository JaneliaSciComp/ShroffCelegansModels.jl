# Copying movies to the OpenShift PVC

Movies (`movie_yz.mp4`, `movie_xz.mp4`) must land at
`/data/annotations/recompute/` on the `shroff-data-test-annotations` PVC
so they are served at `https://shroff-data-test.int.janelia.org/recompute/`.

All running deployment pods mount that PVC **read-only**, so `oc cp` into
any of them is blocked.  The workaround is a temporary busybox pod that
mounts the PVC read-write.

## Step-by-step

```bash
cd deployment

# 1. Spin up a scratch pod with write access to the PVC
pixi run oc run movie-upload --image=busybox --restart=Never \
  --overrides='{
    "spec": {
      "containers": [{
        "name": "movie-upload",
        "image": "busybox",
        "command": ["sh", "-c", "sleep 3600"],
        "volumeMounts": [{"mountPath": "/data", "name": "annotations"}]
      }],
      "volumes": [{
        "name": "annotations",
        "persistentVolumeClaim": {"claimName": "shroff-data-test-annotations"}
      }]
    }
  }' \
  -n shroff-data-test

# 2. Wait for it to be Running
pixi run oc wait pod/movie-upload -n shroff-data-test \
  --for=condition=Ready --timeout=60s

# 3. Copy the movies (adjust local paths as needed)
pixi run oc cp movies/movie_yz.mp4 \
  shroff-data-test/movie-upload:/data/annotations/recompute/movie_yz.mp4
pixi run oc cp movies/movie_xz.mp4 \
  shroff-data-test/movie-upload:/data/annotations/recompute/movie_xz.mp4

# 4. Verify
pixi run oc exec movie-upload -n shroff-data-test -- \
  ls -lh /data/annotations/recompute/movie_yz.mp4 \
         /data/annotations/recompute/movie_xz.mp4

# 5. Delete the scratch pod
pixi run oc delete pod movie-upload -n shroff-data-test
```

## Notes

- The `--overrides` JSON must be a single line or properly escaped in the
  shell; the multi-line form above works with bash here-strings but not
  all shells.
- If the pod image pull is slow, use `oc get pod movie-upload -n shroff-data-test`
  to watch status before copying.
- In production (`shroff-data` namespace) substitute `shroff-data` and
  `shroff-data-annotations` (the production PVC name) accordingly.
- Long-term: `web/scripts/run_recompute_if_needed.jl` already generates
  both movies automatically at the end of each pipeline run, so manual
  copies are only needed when re-running movies outside the normal pipeline.
