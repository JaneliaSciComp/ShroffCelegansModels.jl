using CairoMakie

include(joinpath(@__DIR__, "crop_video.jl"))
include(joinpath(@__DIR__, "meshscatter_average_simple.jl"))

# Headless movie generation using CairoMakie (no display/GPU required).
# Renders all timepoints and re-encodes via crop_video for tight framing.
# avg_dict: Dict from load_latest_average_annotations()
# output_path: destination .mp4 file path
function generate_meshscatter_movie(average_annotations_dict;
        output_path::String, framerate::Int=48, view::Symbol=:yz)
    CairoMakie.activate!()
    raw_path = replace(output_path, ".mp4" => "_raw.mp4")
    with_theme(theme_black()) do
        fig, time_slider, ax = meshscatter_average_simple(average_annotations_dict;
            xy_bounding_radius = -1,
            figure_size = (1920, 360),
            view)
        time_points = time_slider.range[]
        time_slider.blockscene.visible[] = false
        try
            vs = Record(fig, time_points; framerate, px_per_unit=1, preset="fast", update=false) do t
                set_close_to!(time_slider, t)
            end
            save(raw_path, vs)
            crop_video(raw_path, output_path; framerate)
        finally
            isfile(raw_path) && rm(raw_path; force=true)
        end
    end
    bytes = isfile(output_path) ? filesize(output_path) : 0
    @info "Movie written" path=output_path size_mb=round(bytes/1024^2; digits=1) framerate
end
