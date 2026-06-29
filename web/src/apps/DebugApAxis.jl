"""
    DebugApAxis

Submodule for the AP-axis debug viewer (deployment container `debug-ap-axis`,
port 9180). Render via `debug_annotation_ap_axis`; data + render defs come from
the runtime-included `launch_show_average_annotations.jl` (see
[`ShowAverageAnnotations`](@ref) for the pattern). No build-time workload —
shared render compilation is cached by [`CommonScenes`](@ref).
"""
module DebugApAxis

using WGLMakie
using Bonito
using ShroffCelegansModels

module Launch
    using WGLMakie
    using Bonito
    using ShroffCelegansModels
end

const PORT = 9180
const PROXY_PATH = "debug_annotation_ap_axis"

_include_launch() = Base.include(Launch,
    joinpath(pkgdir(ShroffCelegansModels), "scripts", "launch_show_average_annotations.jl"))

function _menu(datasets)
    DOM.div(DOM.ul(map(collect(keys(datasets))) do k
        DOM.li(k), DOM.ul(map(collect(keys(datasets[k]))) do i
            DOM.li(DOM.a(datasets[k][i].path, href="/$PROXY_PATH/$k/$i"))
        end)
    end))
end

function build_server()
    datasets = Launch.datasets
    server = Server("0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/$PROXY_PATH/")
    route!(server, "/" => App(_menu(datasets)))
    for k in keys(datasets)
        for i in keys(datasets[k])
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans debug annotation AP axis") do
                return Launch.debug_annotation_ap_axis(Launch.avg_models, datasets[k][i]; use_myuntwist=true)
            end)
        end
    end
    return server
end

function main()
    @info "Loading data (runtime include)"
    _include_launch()
    WGLMakie.activate!(; resize_to = :body)
    @info "Launching server!" port=PORT
    server = build_server()
    if isinteractive()
        println("Press enter to quit")
        readline()
    else
        wait(server)
    end
end

end # module DebugApAxis
