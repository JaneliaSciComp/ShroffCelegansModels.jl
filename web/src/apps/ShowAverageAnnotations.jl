"""
    ShowAverageAnnotations

Submodule for the average-annotations viewer (deployment container
`show-average-annotations`, port 8180).

The render function `show_average_annotations` and the `datasets`/`avg_models`
globals come from `scripts/launch_show_average_annotations.jl`, which eagerly
reads PVC data at include time. So that include happens at **runtime** inside
[`main`](@ref), into the isolated child module `Launch` (never at precompile
time). Shared Makie/WGLMakie/Bonito render compilation is cached by
[`CommonScenes`](@ref); this submodule carries no build-time workload.
"""
module ShowAverageAnnotations

using WGLMakie
using Bonito
using ShroffCelegansModels

# Child module that will host the runtime-included launch script (render defs +
# datasets/avg_models). Empty at precompile time → no build-time disk reads.
module Launch
    using WGLMakie
    using Bonito
    using ShroffCelegansModels
end

const PORT = 8180
const PROXY_PATH = "show_average_annotations"

_include_launch() = Base.include(Launch,
    joinpath(pkgdir(ShroffCelegansModels), "scripts", "launch_show_average_annotations.jl"))

function build_server()
    datasets = Launch.datasets
    menu = DOM.div(DOM.ul(map(collect(keys(datasets))) do k
        DOM.li(DOM.a(k, href="/$PROXY_PATH/$k"))
    end))
    server = Server("0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/$PROXY_PATH/")
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        route!(server, "/$k" => App(; title="$k: Shroff C. elegans show average annotations") do
            return Launch.show_average_annotations(Launch.avg_models, datasets[k]; use_myuntwist=true)
        end)
    end
    return server
end

function main()
    @info "Loading data (runtime include)"
    _include_launch()
    @info "Priming annotation cache"
    Launch.prime_annotation_caches()
    Launch.alias_cache_unix("/nearline/shroff/")
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

end # module ShowAverageAnnotations
