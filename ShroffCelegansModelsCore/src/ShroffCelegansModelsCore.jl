module ShroffCelegansModelsCore
    using CSV
    using DataFrames
    using GeometryBasics
    using LinearAlgebra
    using BSplineKit
    using QuadGK
    using LRUCache
    using FFTW
    using Dates
    using PrecompileTools: @setup_workload, @compile_workload

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
    # include("parse_worm_dataset_path.jl")

    include("precompile.jl")
end