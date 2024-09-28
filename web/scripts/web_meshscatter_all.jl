using WGLMakie
using Bonito

include("../../src/demo_averaging/save_cache.jl")
include("../../scripts/meshscatter_all.jl")

function meshscatter_all_webapp()
    app = App(; title="Shroff Lab: C. elegans meshscatter_all") do
        return with_theme(meshscatter_all, theme_black())
    end
    server = Server(app, "shroff-data.int.janelia.org", 8082;
        proxy_url="https://shroff-data.int.janelia.org/meshscatter_all/"
    )
    return server
end

function main()
    @info "Loading annotation cache"
    load_annotation_cache()
    @info "Activating WGLMakie"
    WGLMakie.activate!(; resize_to = :body)
    @info "Launching server!"
    meshscatter_all_webapp()
    println("Press enter to quit")
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
