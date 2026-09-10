"""
    FixApAxis

Submodule for the AP-axis fix tool (deployment container `fix-ap-axis`, port
9381). Renders `fix_annotation_ap_axis` with a query-param-driven initial
state and a listener that mirrors annotation/timepoint into the URL. Also runs
the annotation-persist socket server (bound synchronously so bind failures
surface; accept loop spawned).

Data + render globals come from the runtime-included
`launch_show_average_annotations.jl` (see [`ShowAverageAnnotations`](@ref)).
"""
module FixApAxis

using WGLMakie
using Bonito
using Bonito: Session
using HTTP
using URIs
using ShroffCelegansModels

module Launch
    using WGLMakie
    using Bonito
    using ShroffCelegansModels
end

const PORT = 9381
const PROXY_PATH = "fix_annotation_ap_axis"

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
            @info "Creating route /$k/$i for dataset $(datasets[k][i].path)"
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans fix annotation AP axis") do session::Session, request::HTTP.Request
                empty!(ShroffCelegansModels.annotations_cache)
                params = HTTP.URIs.queryparams(URI(request.target).query)
                listener = (a, v) -> evaljs(session, js"""history.replaceState(null, "", "?annotation=" + $a + "&timepoint=$v");""")
                # `fix_annotation_ap_axis` is a `using`-import in the launch
                # script, so it is not reachable as `Launch.fix_annotation_ap_axis`;
                # call it on the package directly (`Launch.avg_models` is a
                # launch-defined binding and stays as-is).
                return ShroffCelegansModels.fix_annotation_ap_axis(
                    Launch.avg_models,
                    datasets[k][i];
                    use_myuntwist=true,
                    initial_timepoint=parse(Int64, get(params, "timepoint", "0")),
                    initial_annotation=get(params, "annotation", nothing),
                    annotation_timepoint_listener=listener,
                )
            end)
        end
    end
    return server
end

function main()
    @info "Loading data (runtime include)"
    _include_launch()
    # The runtime include defines Launch globals at a newer world age than this
    # precompiled `main` can see, so run the rest via `invokelatest`.
    Base.invokelatest() do
        # Bind the persist listener synchronously so a bind failure surfaces here
        # instead of being swallowed by the spawned accept loop.
        persist_listener = ShroffCelegansModels.fix_annotation_ap_axis_persist_listen()
        Threads.@spawn ShroffCelegansModels.fix_annotation_ap_axis_persist_server(persist_listener)
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
end

end # module FixApAxis
