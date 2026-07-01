#!/usr/bin/env julia
#
# Post-deployment liveness check for the Shroff C. elegans Julia web endpoints.
#
# GETs each public endpoint through the nginx front-end and checks two things:
#   1. the HTTP status (a downed upstream shows as a 502/503/504 from nginx), and
#   2. the response body for the signature of a Julia exception — Bonito renders
#      a server-side error into the returned HTML *with HTTP 200*, so a status
#      check alone misses render-time failures (e.g. an UndefVarError thrown
#      while building the page).
#
# Usage:
#   julia --project=web deployment/post_deploy_check.jl [host]
#   julia --project=web deployment/post_deploy_check.jl [host] --poll-until-ready
#
# With --poll-until-ready (alias --time, or SHROFF_CHECK_POLL=1) the script
# starts a clock immediately and polls each endpoint until it first returns
# healthy, reporting per-app time-to-first-healthy — a cold-start measurement.
# Run it right after `oc rollout restart`. Tune with SHROFF_CHECK_POLL_INTERVAL
# (default 2s) and SHROFF_CHECK_POLL_TIMEOUT (default 300s).
#
# Host resolution (first match wins):
#   1. CLI argument                 e.g. shroff-data.int.janelia.org
#   2. SHROFF_CHECK_HOST env var
#   3. SHROFF_HOST env var
#   4. default: shroff-data-test.int.janelia.org
#
# TLS certificate verification is OFF by default (the *.int.janelia.org hosts
# use an internal CA); set SHROFF_CHECK_VERIFY_TLS=1 to turn it on.
#
# Exit code is 0 only if every endpoint is healthy, else 1 — suitable for CI.

using HTTP
using Printf

const DEFAULT_HOST = "shroff-data-test.int.janelia.org"

# (label, path). Root paths keep the trailing slash to match the nginx
# `location /name/` blocks; the dataset page exercises real per-dataset routing
# (which is where a render-time error like the get_datasets_info regression hid).
const ENDPOINTS = [
    ("show_average_annotations",              "/show_average_annotations/"),
    ("show_average_annotations/RW10131",      "/show_average_annotations/RW10131"),
    ("meshscatter_average_edited",            "/meshscatter_average_edited/"),
    ("meshscatter_average_edited_2024_10_24", "/meshscatter_average_edited_2024_10_24/"),
    ("debug_annotation_ap_axis",              "/debug_annotation_ap_axis/"),
    ("debug_annotation_ap_axis_live",         "/debug_annotation_ap_axis_live/"),
    ("debug_annotation_ap_axis_retrack_live", "/debug_annotation_ap_axis_retrack_live/"),
    ("zscore_analysis",                       "/zscore_analysis/"),
    ("modified_times",                        "/modified_times/"),
    ("fix_annotation_ap_axis",                "/fix_annotation_ap_axis/"),
]

# Signatures of a Julia exception / Bonito error page. Kept specific to avoid
# false positives on normal app content.
const ERROR_MARKERS = [
    "Stacktrace", "UndefVarError", "MethodError", "BoundsError", "KeyError",
    "ArgumentError", "DimensionMismatch", "TypeError", "ErrorException",
    "not defined in",
]

# The first non-flag argument is the host; `--`-prefixed args are options.
function target_host()
    pos = filter(a -> !startswith(a, "--"), ARGS)
    !isempty(pos) ? pos[1] :
        get(ENV, "SHROFF_CHECK_HOST", get(ENV, "SHROFF_HOST", DEFAULT_HOST))
end

function error_marker(body)
    for m in ERROR_MARKERS
        occursin(m, body) && return m
    end
    return nothing
end

function probe(url; verify)
    try
        t = @elapsed r = HTTP.get(url; status_exception = false, redirect = true,
                     require_ssl_verification = verify,
                     connect_timeout = 10, readtimeout = 60, retry = false)
        return (r.status, String(r.body), nothing, t)
    catch e
        return (nothing, nothing, e, nothing)
    end
end

is_healthy(status, marker) = status !== nothing && 200 <= status < 400 && marker === nothing

# Single-shot check: probe every endpoint once, PASS/FAIL each. This is the
# steady-state liveness check (run after a deploy has settled).
function single_shot(; base, verify)
    width = maximum(length(first(e)) for e in ENDPOINTS)
    failures = 0
    for (label, path) in ENDPOINTS
        status, body, err, elapsed = probe(base * path; verify)
        marker  = body === nothing ? nothing : error_marker(body)
        healthy = is_healthy(status, marker)
        healthy || (failures += 1)
        mark = healthy ? "PASS" : "FAIL"
        timing = elapsed === nothing ? "  ?.?s" : @sprintf("%5.1fs", elapsed)
        detail = status === nothing ? "ERROR " * sprint(showerror, err) :
                 marker !== nothing ? "HTTP $status, but body contains \"$marker\"" :
                 "HTTP $status"
        println("  $mark  $timing  $(rpad(label, width))  $detail")
    end
    n = length(ENDPOINTS)
    println()
    println(failures == 0 ? "All $n endpoints healthy." : "$failures of $n endpoints FAILED.")
    return failures
end

# Cold-start timing mode: start the clock now (run this immediately after
# `oc rollout restart`), then poll every endpoint until it first returns healthy,
# recording seconds-from-start. This measures per-app container boot + Julia load
# + first-render compile — the cost the PrecompileTools caches are meant to cut.
function poll_until_ready(; base, verify,
        interval = parse(Float64, get(ENV, "SHROFF_CHECK_POLL_INTERVAL", "2")),
        timeout  = parse(Float64, get(ENV, "SHROFF_CHECK_POLL_TIMEOUT", "300")))
    width = maximum(length(first(e)) for e in ENDPOINTS)
    t0 = time()
    ready = Dict{String,Float64}()      # label => seconds-to-first-healthy
    done  = Set{String}()
    println("  clock starts now (run right after `oc rollout restart`); ",
            "interval=$(interval)s timeout=$(timeout)s\n")
    while length(done) < length(ENDPOINTS) && (time() - t0) < timeout
        for (label, path) in ENDPOINTS
            label in done && continue
            status, body, _, _ = probe(base * path; verify)
            marker = body === nothing ? nothing : error_marker(body)
            if is_healthy(status, marker)
                el = time() - t0
                ready[label] = el
                push!(done, label)
                println(@sprintf("  READY  %7.1fs  %s", el, label))
            end
        end
        length(done) < length(ENDPOINTS) && sleep(interval)
    end

    println("\n  Cold-start time-to-first-healthy (seconds from clock start):")
    for (label, _) in ENDPOINTS
        if haskey(ready, label)
            println(@sprintf("    %7.1fs  %s", ready[label], rpad(label, width)))
        else
            println("    TIMEOUT   $(rpad(label, width))  (> $(round(Int, timeout))s)")
        end
    end
    timed_out = length(ENDPOINTS) - length(done)
    println()
    println(timed_out == 0 ?
        "All $(length(ENDPOINTS)) endpoints became healthy." :
        "$timed_out of $(length(ENDPOINTS)) endpoints never became healthy within $(round(Int, timeout))s.")
    return timed_out
end

function main()
    host   = target_host()
    verify = get(ENV, "SHROFF_CHECK_VERIFY_TLS", "0") == "1"
    base   = "https://$host"
    poll   = any(a -> a in ("--poll-until-ready", "--time"), ARGS) ||
             get(ENV, "SHROFF_CHECK_POLL", "0") == "1"

    println("Post-deployment endpoint check$(poll ? " — cold-start timing" : "")")
    println("  base: $base")
    println("  TLS verification: $(verify ? "on" : "off")\n")

    failures = poll ? poll_until_ready(; base, verify) : single_shot(; base, verify)
    exit(failures == 0 ? 0 : 1)
end

main()
