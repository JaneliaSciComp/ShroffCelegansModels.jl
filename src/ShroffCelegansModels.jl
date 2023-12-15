module ShroffCelegansModels
	using CSV
	#using Makie
	using Makie
	using DataFrames
	using GeometryBasics
	#using ColorSchemes
	using LinearAlgebra
	using BSplineKit
	using QuadGK
	using LRUCache
	using PrecompileTools: @setup_workload, @compile_workload

	include("util.jl")
	include("ParametricSplines.jl")
	include("types.jl")

	using .ParametricSplines
	using .Types

	include("build_model.jl")
	include("mesh.jl")
	include("straighten.jl")

	include("average.jl")
	include("datasets.jl")

	using .Datasets

	include("show.jl")
	include("points.jl")

	using .Points

	include("model_time_series.jl")

	include("precompile.jl")

end
