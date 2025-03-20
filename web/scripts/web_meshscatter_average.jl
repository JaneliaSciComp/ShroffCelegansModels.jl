using WGLMakie
using Bonito
using ShroffCelegansModels

#push!(LOAD_PATH, "/groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl")
push!(LOAD_PATH, dirname(dirname(pathof(ShroffCelegansModels))))
include("../../src/demo_averaging/save_cache.jl")
include("../../src/demo_averaging/average_annotations.jl")
include("../../scripts/meshscatter_average.jl")

function black_body(fig)
    DOM.body(fig, style=Styles(CSS("background-color" => "black")))
end

function meshscatter_average_webapp()
    average_annotation_dict = load_average_annotations()
    app = App(; title="Shroff Lab: C. elegans meshscatter_average") do
        return with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict))
        end
    end
    nerve_ring_app = App(; title="Shroff Lab: C. elegans meshscatter_average/nerve_ring") do
        return with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict; nerve_ring=true))
        end
    end
    server = Server(app, "shroff-data.int.janelia.org", 8380;
        proxy_url="https://shroff-data.int.janelia.org/meshscatter_average/"
    )
    route!(server, "/nerve_ring" => nerve_ring_app)
    return server
end

function main()
    @info "Loading annotation cache"
    load_annotation_cache()
    @info "Activating WGLMakie"
    WGLMakie.activate!(; resize_to = :body)
    @info "Launching server!"
    meshscatter_average_webapp()
    println("Press enter to quit")
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
