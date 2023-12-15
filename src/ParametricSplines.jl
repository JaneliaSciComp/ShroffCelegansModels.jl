module ParametricSplines
    using BSplineKit: BSplineKit, knots, Derivative, BSplineOrder, Natural, SplineInterpolation
    using QuadGK: QuadGK, quadgk
    using DataFrames: DataFrame
    using GeometryBasics: Point, Point3
    using LinearAlgebra: norm

    #export NDParametricSplineInterpolation
    export interpolate_natural_cubic_spline

	#=
	struct NDParametricSplineInterpolation{T}
		splines::Vector{T}
	end
	function (itps::NDParametricSplineInterpolation)(time)
		Point{length(itps.splines)}(map(itps.splines) do s
			s(time)
		end)
	end
	function Base.getindex(itps::NDParametricSplineInterpolation, elements...)
		getindex(itps.splines, elements...)
	end
	function Base.:*(d::Derivative, s::NDParametricSplineInterpolation)
		NDParametricSplineInterpolation(d .* s.splines)
	end
	BSplineKit.knots(s::NDParametricSplineInterpolation) = knots(first(s.splines))
	function QuadGK.quadgk(s::NDParametricSplineInterpolation)
		d = Derivative(1) * s
		quadgk(unique(knots(s))...) do t
			norm(d(t))
		end
	end
	function QuadGK.quadgk()
	end
	=#

	function QuadGK.quadgk(s::SplineInterpolation)
		d = Derivative(1) * s
		quadgk(unique(knots(s))...) do t
			norm(d(t))
		end
	end
# ╔═╡ d323b169-b352-46ad-9d35-a4da784cc710
"""
	interpolate_natural_cubic_spline(time, data)
"""
function interpolate_natural_cubic_spline(time, data::Union{AbstractMatrix, DataFrame})
	#=
	itp_columns = interpolate_natural_cubic_spline.((time,), eachcol(data))
	NDParametricSplineInterpolation(itp_columns)
	=#
	return interpolate_natural_cubic_spline(time, Point3.(Vector.(eachrow(data))))
end

#=
# ╔═╡ 96d3a2f4-93ce-43aa-99d8-2fc1fbb21b1e
function interpolate_natural_cubic_spline(time, data::Vector{<: Point3})
	interpolate_natural_cubic_spline(time, stack(data; dims=1))
end
=#

# ╔═╡ 86b367d5-506d-4ed8-a3ca-a356e4b5be1a
function interpolate_natural_cubic_spline(time, data::AbstractVector)
	return BSplineKit.interpolate(time, data, BSplineOrder(4), Natural())
end    

end