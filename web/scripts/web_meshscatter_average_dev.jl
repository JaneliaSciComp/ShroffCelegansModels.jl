using WGLMakie
using Bonito
using ShroffCelegansModels
using Sockets

include("../../scripts/meshscatter_average_dev.jl")

using ShroffCelegansModels: load_average_annotations, load_latest_average_annotations, load_annotation_cache

function black_body(fig)
    DOM.body(fig, style=Styles(CSS("background-color" => "black")))
end

function meshscatter_average_webapp()
    average_annotation_dict = load_latest_average_annotations(
        default_filename = "edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells_2025_02_13.h5",
    )
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

    shroff_data_ip = "0.0.0.0"
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
