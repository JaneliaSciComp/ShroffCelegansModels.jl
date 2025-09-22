module Types
    using GeometryBasics: Point, Point3f, Point3
	using QuadGK: quadgk
	using LinearAlgebra: norm
	using BSplineKit: Derivative, knots
	import BSplineKit.SplineInterpolations: interpolation_points

	using ..ParametricSplines: interpolate_natural_cubic_spline

	export AbstractCelegansModel, CelegansModel, CachedCelegansModel, StraightenedCelegansModel
	export central_spline, central_spline_d1
	export transverse_spline, num_transverse_splines, transverse_splines
	export lattice_knots
	export model_length
	export names
	export interpolation_points

    abstract type AbstractCelegansModel{S} end
    central_spline_d1(model::AbstractCelegansModel) = Derivative(1) * central_spline(model)
    model_length(model::AbstractCelegansModel) = first(quadgk(central_spline(model)))
    length(model::AbstractCelegansModel) = ceil(Int, model_length(model)) + 1
    Base.length(model::AbstractCelegansModel) = length(model)
	transverse_splines(model::AbstractCelegansModel) = map(1:num_transverse_splines(model)) do i
		transverse_spline(model, i)
	end
	lattice_knots(model::AbstractCelegansModel) = knots(central_spline(model))[4:end-3]
	interpolation_points(model::AbstractCelegansModel) = interpolation_points(first(transverse_splines(model)))

    struct CelegansModel{S, N} <: AbstractCelegansModel{S}
        transverse_splines::Vector{S}
        central_spline::S
		names::Vector{N}
    end

    central_coordinate(model::CelegansModel, t::Float64) = model.central_spline(t)
	central_spline(model::CelegansModel) = model.central_spline
	transverse_spline(model::CelegansModel, i::Integer) = model.transverse_splines[i]
	num_transverse_splines(model::CelegansModel) = Base.length(model.transverse_splines)
	names(model::CelegansModel) = model.names

	struct CachedCelegansModel{S} <: AbstractCelegansModel{S}
		model::CelegansModel{S}
		central_spline_d1::S
		model_length::Float64
		length::Int
		function CachedCelegansModel(model::AbstractCelegansModel{S})  where S
			new{S}(
				model,
				central_spline_d1(model),
				model_length(model),
				length(model)
			)
		end
	end
	central_coordinate(model::CachedCelegansModel, t::Float64) = central_coordinate(model.model, t)
	central_spline(model::CachedCelegansModel) = central_spline(model.model)
	transverse_spline(model::CachedCelegansModel, i::Integer) = transverse_spline(model.model, i)
	central_spline_d1(model::CachedCelegansModel) = model.central_spline_d1
	model_length(model::CachedCelegansModel) = model.model_length
	length(model::CachedCelegansModel) = model.length
    Base.length(model::CachedCelegansModel) = length(model)
	parent(model::CachedCelegansModel) = model.model
	Base.parent(model::CachedCelegansModel) = parent(model)
	num_transverse_splines(model::CachedCelegansModel) = num_transverse_splines(parent(model))
	names(model::CachedCelegansModel) = names(parent(model))



	struct StraightenedTransverseSpline{S}
		transverse_spline::S
		central_spline::S
		straightened_central_spline::S
		transverse_vector::Point3f
	end
	function (s::StraightenedTransverseSpline)(t)
		r = norm(s.transverse_spline(t) - s.central_spline(t))
		return s.straightened_central_spline(t) + r*s.transverse_vector
	end
	function StraightenedTransverseSpline(model, straightened_model, i::Integer)
		degree = 360/Base.length(model.transverse_splines)*(i-1)
		StraightenedTransverseSpline(
			model.transverse_splines[i],
			model.central_spline,
			straightened_model.central_spline,
			Point3f(cosd(degree), sind(degree), 0)
		)
	end
	interpolation_points(sts::StraightenedTransverseSpline) = interpolation_points(sts.transverse_spline)

	struct StraightenedCelegansModel{S} <: AbstractCelegansModel{S}
		twisted_model::CelegansModel{S}
		central_spline::S
	end
	function StraightenedCelegansModel(twisted_model::CelegansModel{S}) where S
		d = central_spline_d1(twisted_model)
		r = LinRange(0, 1, length(twisted_model))
		n = norm.(d.(r))
		s = cumsum([0; (n[1:end-1] .+ n[2:end])./2]) .* step(r)
		#z = interpolate_natural_cubic_spline([0,0.25,0.5,0.75,1],[0.0,0.0,0.0,0.0,0.0])
		s = map(s) do z
			Point3(0, 0, z)
		end
		central_straightened = interpolate_natural_cubic_spline(collect(r), s)
		#S([z, z, central_straightened])
		#typeof(central_straighted)
		StraightenedCelegansModel(
			twisted_model,
			central_straightened
		)
	end
	function transverse_spline(smodel::StraightenedCelegansModel, i)
		return StraightenedTransverseSpline(smodel.twisted_model, smodel, i)
	end
	parent(smodel::StraightenedCelegansModel) = smodel.twisted_model
	Base.parent(smodel::StraightenedCelegansModel) = parent(smodel)
	central_spline(smodel::StraightenedCelegansModel) = smodel.central_spline
	length(smodel::StraightenedCelegansModel) = length(parent(smodel))
	num_transverse_splines(model::StraightenedCelegansModel) = num_transverse_splines(parent(model))
	names(model::StraightenedCelegansModel) = names(parent(model))
	lattice_knots(model::StraightenedCelegansModel) = lattice_knots(parent(model))

	"""
		RadialCelegansModel

	The radial C. elegans model consists of a 3-dimensional parametric spline describing the central axis of the worm.
	The transverse axes that run from anterior to posterior are defined by their radial distance from the the central spline
	along regular radial spokes normal to the derivative of the central spline. These are modeled as 1-dimensional interpolate_natural_cubic_spline
	cubic splines.

	This differs from `CelegansModel` where the transverse_splines are 3-dimensional parameteric splines rather than 1-dimensional splines.
	"""
	struct RadialCelegansModel{C,R} <: AbstractCelegansModel{C}
		transverse_splines::Vector{R}
		central_spline::C
		names::Vector{String}
	end

end