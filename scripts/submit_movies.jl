#!/usr/bin/env julia
# Submit GLMakie GPU movie jobs to LSF gpu_l4.
#
# Usage:
#   julia scripts/submit_movies.jl
#   julia scripts/submit_movies.jl /path/to/avg.h5
#   julia scripts/submit_movies.jl /path/to/avg.h5 /path/to/output_dir
#   julia scripts/submit_movies.jl /path/to/avg.h5 /path/to/output_dir yz xz xy
#
# Default views: yz xz
# Output files: <output_dir>/movie_<view>.mp4

const REPO   = abspath(joinpath(@__DIR__, ".."))
const JULIA  = joinpath(homedir(), ".juliaup", "bin", "julia")
const SCRIPT = joinpath(REPO, "glmakie", "scripts", "generate_movie_glmakie.jl")

const H5_DEFAULT = joinpath(REPO, "movies",
    "edited_smoothed_average_annotations_r020_theta020_z030_2026_06_09_110120.h5")

h5_path    = get(ARGS, 1, H5_DEFAULT)
output_dir = get(ARGS, 2, joinpath(REPO, "movies"))
views      = length(ARGS) >= 3 ? ARGS[3:end] : ["yz", "xz"]

mkpath(output_dir)

for view in views
    out     = joinpath(output_dir, "movie_$(view).mp4")
    log     = joinpath(output_dir, "movie_$(view)_%J.log")
    err_log = joinpath(output_dir, "movie_$(view)_%J.err")
    env_str = "all,H5=$(h5_path),OUT=$(out),VIEW=$(view)"
    gpu_str = "num=1:j_exclusive=yes"
    mem_str = "rusage[mem=32000]"
    job_cmd = "time xvfb-run -a $(JULIA) --project=$(REPO)/glmakie $(SCRIPT) \$H5 \$OUT \$VIEW"

    cmd = `bsub
        -J shroff_movie_$(view)
        -P scicompsoft
        -q gpu_l4
        -gpu $(gpu_str)
        -n 4
        -R $(mem_str)
        -W 00:30
        -o $(log)
        -e $(err_log)
        -env $(env_str)
        $(job_cmd)`

    @info "Submitting $(view) view" output=out h5=h5_path
    run(cmd)
end
