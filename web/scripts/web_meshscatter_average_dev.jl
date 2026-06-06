using WGLMakie
using Bonito
using ShroffCelegansModels
using Sockets

#push!(LOAD_PATH, "/groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl")
push!(LOAD_PATH, dirname(dirname(pathof(ShroffCelegansModels))))
include("../../src/demo_averaging/save_cache.jl")
include("../../src/demo_averaging/average_annotations.jl")
include("../../scripts/meshscatter_average_dev.jl")

function black_body(fig)
    DOM.body(fig, style=Styles(CSS("background-color" => "black")))
end

function meshscatter_average_webapp()
    average_annotation_dict = load_average_annotations()
    app = App(; title="Shroff Lab: C. elegans meshscatter_average") do session::Session
        return with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict; session))
        end
    end
    nerve_ring_app = App(; title="Shroff Lab: C. elegans meshscatter_average/nerve_ring") do
        return with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict; nerve_ring=true))
        end
    end

    shroff_data_ip = "0.0.0.0",
    server = Server(app, shroff_data_ip, 8480;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/meshscatter_average_dev/"
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
    server = meshscatter_average_webapp()
    if isinteractive()
        println("Press enter to quit")
        readline()
    else
        wait(server)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
