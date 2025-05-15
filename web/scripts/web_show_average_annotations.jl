using WGLMakie
using Bonito


if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

includet("../../scripts/launch_show_average_annotations.jl")

function web_show_average_annotations(datasets = datasets)
    menu = DOM.div(
        DOM.ul(
            map(collect(keys(datasets))) do k
                DOM.li(DOM.a(k, href="/show_average_annotations/$k"))
            end
        )
    )
    server = Server(
        "shroff-data.int.janelia.org", 8180;
        proxy_url="https://shroff-data.int.janelia.org/show_average_annotations/"
    )
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        route!(server, "/$k" => App(; title="$k: Shroff C. elegans show average annotations") do
                return show_average_annotations(avg_models, datasets[k]; use_myuntwist=true);
        end)
    end
    return server
end

function web_main()
    WGLMakie.activate!(; resize_to = :body)
    server = web_show_average_annotations()
end

println(@__FILE__)
println(abspath(PROGRAM_FILE))

if run_web_main
    wait(web_main())
    println("Press enter to quit")
    readline()
end
