# Cold-start benchmark for the meshscatter_average web app.
#
# Run in a FRESH Julia process so the timings reflect what a container sees on
# startup:
#
#     julia --project=web web/scripts/bench_cold_start.jl
#
# With the PrecompileTools workload baked into the package image, the first
# `build_apps` + `export_static` call should be a small fraction of what it
# costs without the cache. Compare by toggling the workload and reprecompiling.
using ShroffCelegansModelsWebInterface: MeshscatterAverage
using Bonito: export_static, App, Session
using WGLMakie

const M = MeshscatterAverage

function main()
    dict = M._synthetic_annotation_dict()

    t_load = @elapsed WGLMakie.activate!()
    @info "WGLMakie.activate!" seconds = t_load

    t_build = @elapsed app, _ = M.build_apps(dict)
    @info "first build_apps()" seconds = t_build

    t_render = @elapsed mktempdir() do dir
        export_static(joinpath(dir, "out.html"), app)
    end
    @info "first export_static()" seconds = t_render

    @info "TOTAL time-to-first-render" seconds = t_load + t_build + t_render
end

main()
