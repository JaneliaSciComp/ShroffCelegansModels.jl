using GLMakie
using ShroffCelegansModels

const _repo_root = joinpath(@__DIR__, "..", "..")
const _scripts_dir = joinpath(_repo_root, "scripts")

include(joinpath(_scripts_dir, "crop_video.jl"))
include(joinpath(_scripts_dir, "meshscatter_average_simple.jl"))

# GLMakie offscreen rendering via EGL (no X display required).
# Set JULIA_GLMAKIE_BACKEND=egl in the environment before loading, or
# export the variable in the bsub script.  When a real X display is
# available (e.g. xvfb-run) the default OpenGL backend is used automatically.
function generate_meshscatter_movie_gl(average_annotations_dict;
        output_path::String,
        framerate::Int = 48,
        view::Symbol = :yz)
    raw_path = replace(output_path, ".mp4" => "_raw.mp4")
    with_theme(theme_black()) do
        fig, time_slider, ax = meshscatter_average_simple(average_annotations_dict;
            xy_bounding_radius = -1,
            figure_size = (1920, 360),
            show_legend = true,
            view = view)
        time_points = time_slider.range[]
        time_slider.blockscene.visible[] = false
        try
            # Record with update=false preserves the zoom set above;
            # the default update=true calls reset_limits! which resets the camera.
            vs = Record(fig, time_points; framerate, px_per_unit = 1, update = false) do t
                set_close_to!(time_slider, t)
            end
            save(raw_path, vs)
            crop_video(raw_path, output_path; framerate)
        finally
            isfile(raw_path) && rm(raw_path; force = true)
        end
    end
    bytes = isfile(output_path) ? filesize(output_path) : 0
    @info "Movie written" path=output_path size_mb=round(bytes/1024^2; digits=1) framerate
end

# ── CLI entry point ────────────────────────────────────────────────────────────
# Usage: julia script.jl [h5_path [output_path [view]]]
# view: yz (default), xz, or xy
h5_path = length(ARGS) >= 1 ? ARGS[1] :
    "/tmp/edited_smoothed_average_annotations_r020_theta020_z030_2026_06_09_110120.h5"

output_path = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(_repo_root, "movies", "test_movie_glmakie.mp4")

view_arg = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :yz

mkpath(dirname(output_path))

@info "Loading average annotations" h5_path
avg_dict = ShroffCelegansModels.load_latest_average_annotations(
    default_filename = basename(h5_path),
    dir = dirname(h5_path),
)
@info "Loaded" n_groups = length(avg_dict)

GLMakie.activate!()
@info "Generating movie (GLMakie)" output_path view=view_arg
generate_meshscatter_movie_gl(avg_dict; output_path, view = view_arg)
