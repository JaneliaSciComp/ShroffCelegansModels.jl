"""
Dispatch a one-shot recompute Job to the in-cluster Kubernetes API.

The Job spec is built from environment variables so the spec can change without
rebuilding the image:

  RECOMPUTE_JOB_IMAGE      — container image (defaults to env IMAGE, then a placeholder)
  RECOMPUTE_JOB_COMMAND    — JSON-encoded array, e.g. ["julia","--version"]
  RECOMPUTE_JOB_NAME_PREFIX — Job name prefix (default "recompute-averages")
  RECOMPUTE_JOB_BACKOFF_LIMIT, RECOMPUTE_JOB_ACTIVE_DEADLINE_SECONDS — passthroughs
  RECOMPUTE_JOB_TTL_SECONDS — spec.ttlSecondsAfterFinished (default 86400)

Auth uses the pod's projected ServiceAccount token at
/var/run/secrets/kubernetes.io/serviceaccount/. If that path is absent (running
outside K8s, e.g. local dev) the call is a no-op that logs a warning.
"""

using HTTP: HTTP
using ShroffCelegansModels.JSON3
using Dates: now, format

const SA_DIR = "/var/run/secrets/kubernetes.io/serviceaccount"

# TLS verification for https://kubernetes.default.svc is intentionally scoped
# off at the CronJob env level via JULIA_SSL_NO_VERIFY_HOSTS=kubernetes.default.svc
# (HTTP.jl reads this via NetworkOptions.verify_host). The bearer token is what
# authenticates the request; the in-cluster API endpoint is reachable only
# through cluster networking, so cert verification adds little here while
# requiring us to build a custom SSL context per call.

function _read_sa_file(name::AbstractString)::Union{Nothing, String}
    path = joinpath(SA_DIR, name)
    isfile(path) || return nothing
    return strip(read(path, String))
end

function _build_job_spec(timestamp::AbstractString, namespace::AbstractString, change_count::Int)
    image = get(ENV, "RECOMPUTE_JOB_IMAGE") do
        get(ENV, "IMAGE", "image-registry.openshift-image-registry.svc:5000/$(namespace)/$(namespace):latest")
    end
    command_json = get(ENV, "RECOMPUTE_JOB_COMMAND", "[\"julia\",\"--version\"]")
    command = JSON3.read(command_json, Vector{String})
    name_prefix = get(ENV, "RECOMPUTE_JOB_NAME_PREFIX", "recompute-averages")
    backoff_limit = parse(Int, get(ENV, "RECOMPUTE_JOB_BACKOFF_LIMIT", "1"))
    active_deadline = parse(Int, get(ENV, "RECOMPUTE_JOB_ACTIVE_DEADLINE_SECONDS", "14400"))
    ttl_seconds = parse(Int, get(ENV, "RECOMPUTE_JOB_TTL_SECONDS", "86400"))

    return Dict(
        "apiVersion" => "batch/v1",
        "kind" => "Job",
        "metadata" => Dict(
            "generateName" => "$(name_prefix)-",
            "namespace" => namespace,
            "labels" => Dict(
                "app" => namespace,
                "component" => name_prefix,
                "triggered-by" => "save-modified-times",
            ),
            "annotations" => Dict(
                "shroff-data/triggered-at" => timestamp,
                "shroff-data/mtime-changes" => string(change_count),
            ),
        ),
        "spec" => Dict(
            "backoffLimit" => backoff_limit,
            "activeDeadlineSeconds" => active_deadline,
            "ttlSecondsAfterFinished" => ttl_seconds,
            "template" => Dict(
                "metadata" => Dict("labels" => Dict(
                    "app" => namespace,
                    "component" => name_prefix,
                )),
                "spec" => Dict(
                    "restartPolicy" => "Never",
                    "containers" => [Dict(
                        "name" => name_prefix,
                        "image" => image,
                        "command" => command,
                        "workingDir" => "/app/ShroffCelegansModels.jl",
                    )],
                ),
            ),
        ),
    )
end

"""
    dispatch_recompute_job(change_count) -> Union{Nothing, String}

POST a recompute Job to the in-cluster K8s API. Returns the created Job name on
success, `nothing` if dispatch was skipped (no SA token / disabled). Raises if
the API call itself fails — the caller decides whether to swallow.
"""
function dispatch_recompute_job(change_count::Int)::Union{Nothing, String}
    token = _read_sa_file("token")
    namespace = _read_sa_file("namespace")
    if token === nothing || namespace === nothing
        @warn "Not running in a Kubernetes pod (no ServiceAccount token) — skipping dispatch" SA_DIR
        return nothing
    end

    timestamp = format(now(), "yyyy-mm-ddTHH:MM:SS")
    spec = _build_job_spec(timestamp, namespace, change_count)
    body = JSON3.write(spec)

    url = "https://kubernetes.default.svc/apis/batch/v1/namespaces/$(namespace)/jobs"
    headers = [
        "Authorization" => "Bearer $(token)",
        "Content-Type" => "application/json",
        "Accept" => "application/json",
    ]

    @info "Dispatching recompute Job" url namespace change_count
    resp = HTTP.post(url, headers, body; status_exception=true)
    parsed = JSON3.read(resp.body)
    job_name = get(get(parsed, :metadata, Dict()), :name, "<unknown>")
    @info "Recompute Job created" job_name
    return String(job_name)
end
