using WGLMakie
using Bonito
using Revise


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
        string(Sockets.getaddrinfo("shroff-data.int.janelia.org")), 9381;
        proxy_url="https://shroff-data.int.janelia.org/fix_annotation_ap_axis/"
    )
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        for i in keys(datasets[k])
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans fix annotation AP axis") do
                empty!(annotations_cache)
                return fix_annotation_ap_axis(avg_models, datasets[k][i]; use_myuntwist=true);
            end)
        end
    end
    return server
end

function web_main()
    Threads.@spawn fix_annotation_ap_axis_persist_server()
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
