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

    include("precompile.jl")

end
