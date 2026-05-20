module ShroffCelegansModels
    using BSplineKit
    using CSV
    using ColorSchemes: ColorSchemes
    using Colors: Colors
    using CoordinateTransformations: CoordinateTransformations
    using DataFrames
    using Dates
    using FFTW
    using FileIO: FileIO
    using FixedPointNumbers: FixedPointNumbers
    using GeometryBasics
    using HDF5: HDF5
    using Interpolations: Interpolations
    using JSON3: JSON3
    using LRUCache
    using LinearAlgebra
    using Makie
    using Missings: Missings
    using Observables: Observables
    using Pkg: Pkg
    using PrecompileTools: @setup_workload, @compile_workload
    using ProgressMeter: ProgressMeter
    using QuadGK
    using Sockets: Sockets
    using StaticArrays: StaticArrays
    using Statistics: Statistics
    using StatsBase: StatsBase
    using ThinPlateSplines: ThinPlateSplines
    using TiffImages: TiffImages
    if gethostname() == "KITTISOPIKULM-2"
        const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_2026_03_19_v2.json"
    else
        const config_path = joinpath(@__DIR__, "..", "config", "linux", "config_2026_03_19_v2.json")
    end
    const voxel_size = 0.1625 # um

    include("datasets.jl")
    include("MIPAVIO.jl")   

    include("util.jl")
    include("ParametricSplines.jl")
    include("types.jl")

    using .ParametricSplines
    using .Types

    include("build_model.jl")
    include("mesh.jl")
    include("straighten.jl")

    include("average.jl")

    using .Datasets

    include("show.jl")
    include("points.jl")

    using .Points

    include("model_time_series.jl")
    include("area.jl")

    include("annotation_untwist.jl")
    include("parse_worm_dataset_path.jl")

    include("demo_averaging/read_config_json.jl")
    # save_celegans_avg_models loads modelio
    include("save_celegans_avg_models.jl")
    include("demo_averaging/seam_cell_pts.jl")
    include("demo_averaging/fix_annotation_ap_axis.jl")
    include("demo_averaging/get_avg_models.jl")
    include("demo_averaging/transform_annotations.jl")
    include("demo_averaging/get_group_annotation_positions_over_time.jl")
    include("demo_averaging/average_annotations.jl")
    include("makie.jl")

    include("precompile.jl")

end
