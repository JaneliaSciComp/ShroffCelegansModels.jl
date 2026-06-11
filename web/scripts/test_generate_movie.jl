using ShroffCelegansModels

const _scripts_dir = joinpath(@__DIR__, "..", "..", "scripts")
include(joinpath(_scripts_dir, "generate_pipeline_movie.jl"))

h5_path = length(ARGS) >= 1 ? ARGS[1] :
    "/tmp/edited_smoothed_average_annotations_r020_theta020_z030_2026_06_09_110120.h5"

output_path = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "..", "..", "movies", "test_movie.mp4")

mkpath(dirname(output_path))

@info "Loading average annotations" h5_path
avg_dict = ShroffCelegansModels.load_latest_average_annotations(
    default_filename = basename(h5_path),
    dir = dirname(h5_path),
)
@info "Loaded" n_groups=length(avg_dict)

@info "Generating movie" output_path
generate_meshscatter_movie(avg_dict; output_path)
