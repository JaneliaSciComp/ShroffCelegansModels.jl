using WGLMakie
using Bonito
using Revise


if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

includet("../../scripts/launch_show_average_annotations.jl")
#includet("../../src/demo_averaging/show_average_annotations.jl")

function web_show_average_annotations(datasets = datasets)
    menu = DOM.div(
        DOM.ul(
            map(collect(keys(datasets))) do k
                DOM.li(DOM.a(k, href="/show_average_annotations/$k"))
            end
        )
    )
    shroff_data_ip = "0.0.0.0"
    server = Server(
        shroff_data_ip, 8180;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/show_average_annotations/"
    )
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        route!(server, "/$k" => App(; title="$k: Shroff C. elegans show average annotations") do
                return show_average_annotations(avg_models, datasets[k]; use_myuntwist=true);
                #Revise.retry()
                #return @invokelatest show_average_annotations(avg_models, datasets[k]; use_myuntwist=true);
        end)
    end
    return server
end

function web_main()
    # Prime both annotation caches at runtime so the first render reuses
    # precomputed positions instead of recomputing (group positions via
    # my_annotation_position_cache; untwisted annotations via annotations_cache).
    # This is done explicitly here because the module no longer loads them at
    # precompile (that froze the build-time snapshot into the .ji). alias_cache_unix
    # then maps the Windows-rooted keys to the Linux dataset.path used at runtime.
    @info "Loading annotation cache"
    prime_annotation_caches()
    alias_cache_unix("/nearline/shroff/")
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
