using WGLMakie
using Bonito
using Revise
using HTTP
using URIs

if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

include("../../scripts/launch_show_average_annotations.jl")

function web_debug_annotation_ap_axis(datasets = datasets)
    menu = DOM.div(
        DOM.ul(
            map(collect(keys(datasets))) do k
                DOM.li(k),
                DOM.ul(
                    map(collect(keys(datasets[k]))) do i
                        DOM.li(DOM.a(datasets[k][i].path, href="/fix_annotation_ap_axis/$k/$i"))
                    end
                )
            end
        )
    )
    server = Server(
        "0.0.0.0", 9381;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/fix_annotation_ap_axis/"
    )
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        for i in keys(datasets[k])
            @info "Creating routing /$k/$i for dataset $(datasets[k][i].path)"
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans fix annotation AP axis") do session::Bonito.Session, request::HTTP.Request
                empty!(annotations_cache)
                params = HTTP.URIs.queryparams(URI(request.target).query)
                listener = (a,v) -> evaljs(session, js"history.replaceState(null, \"\", \"?annotation=\" + $a + \"&timepoint=$v\");")
                return fix_annotation_ap_axis(
                    avg_models,
                    datasets[k][i];
                    use_myuntwist=true,
                    initial_timepoint=parse(Int64, get(params, "timepoint", "0")),
                    initial_annotation=get(params, "annotation", nothing),
                    annotation_timepoint_listener=listener
                );
            end)
        end
    end
    return server
end

function web_main()
    # Bind the persist listener synchronously so a bind failure surfaces here
    # instead of being silently swallowed by the spawned task (see
    # fix_annotation_ap_axis_persist_listen). Only the accept loop is spawned.
    persist_listener = fix_annotation_ap_axis_persist_listen()
    Threads.@spawn fix_annotation_ap_axis_persist_server(persist_listener)
    WGLMakie.activate!(; resize_to = :body)
    server = web_debug_annotation_ap_axis()
end

println(@__FILE__)
println(abspath(PROGRAM_FILE))

if run_web_main
    wait(web_main())
    println("Press enter to quit")
    readline()
end
