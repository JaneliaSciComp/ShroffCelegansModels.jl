using WGLMakie
using Bonito

include(joinpath(@__DIR__, "meshscatter_average_simple.jl"))

# Export a self-contained interactive HTML file via Bonito.export_static.
# The resulting HTML embeds all WGLMakie JavaScript and cell-position data so
# the time slider works in any browser without a Julia backend.
# avg_dict: Dict from load_latest_average_annotations()
# output_path: destination .html file path
function export_meshscatter_static(average_annotations_dict; output_path::String)
    WGLMakie.activate!()
    app = App(; title="C. elegans Average Annotations") do session::Session
        with_theme(theme_black()) do
            fig, _ = meshscatter_average_simple(average_annotations_dict; xy_bounding_radius=sqrt(52))
            DOM.body(fig, style=Styles(CSS("background-color" => "black")))
        end
    end
    Bonito.export_static(output_path, app)
    bytes = isfile(output_path) ? filesize(output_path) : 0
    @info "Static export written" path=output_path size_mb=round(bytes/1024^2; digits=1)
end
