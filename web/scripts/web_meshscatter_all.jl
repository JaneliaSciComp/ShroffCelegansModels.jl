using WGLMakie
using Bonito
using ShroffCelegansModels

include("../../scripts/meshscatter_all.jl")

using ShroffCelegansModels: load_annotation_cache

function meshscatter_all_webapp()
    app = App(; title="Shroff Lab: C. elegans meshscatter_all") do
        return with_theme(meshscatter_all, theme_black())
    end
    server = Server(app, "0.0.0.0", 8082;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/meshscatter_all/"
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
