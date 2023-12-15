# ╔═╡ 9d9e7f86-dc98-4d5b-a0cd-e6a19af34e45
function straighten_celegans_model(model::CelegansModel)
	r = LinRange(0, 1, length(model))
	s = let central_spline_d1 = Derivative(1) * model.central_spline
		dr = diff(r)
		s = cumsum(norm.(central_spline_d1.(r)) .* dr[1])
	end
	sections_radii = map(r) do x
		map(model.transverse_splines) do s
			norm(s(x) - model.central_spline(x))
		end
	end
	normal_vector = Point3f(0.0, 0.0, 1.0)
	right_vector = Point3f(1.0, 0.0, 0.0)
	ap_vector = Point3f(0.0, 1.0, 0.0)
	circ_pts = get_circle_points(right_vector, normal_vector, Point3f(0.0, 0.0, 0.0))
	map(zip(sections_radii,s)) do (radii, ds)
		Point3f.(circ_pts .* radii .+ (ap_vector*ds,))
	end
end