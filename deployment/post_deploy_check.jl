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

target_host() = !isempty(ARGS) ? ARGS[1] :
    get(ENV, "SHROFF_CHECK_HOST", get(ENV, "SHROFF_HOST", DEFAULT_HOST))

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

function main()
    host   = target_host()
    verify = get(ENV, "SHROFF_CHECK_VERIFY_TLS", "0") == "1"
    base   = "https://$host"
    width  = maximum(length(first(e)) for e in ENDPOINTS)

    println("Post-deployment endpoint check")
    println("  base: $base")
    println("  TLS verification: $(verify ? "on" : "off")\n")

    failures = 0
    for (label, path) in ENDPOINTS
        url = base * path
        status, body, err, elapsed = probe(url; verify)
        marker  = body === nothing ? nothing : error_marker(body)
        healthy = status !== nothing && 200 <= status < 400 && marker === nothing
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
    if failures == 0
        println("All $n endpoints healthy.")
    else
        println("$failures of $n endpoints FAILED.")
    end
    exit(failures == 0 ? 0 : 1)
end

main()
