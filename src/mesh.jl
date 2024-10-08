	function get_model_contour_mesh(model_contours_path::AbstractString)
		sections = get_model_contour_sections(model_contours_path)
		get_model_contour_mesh(sections)
	end
	function get_model_contour_mesh(
		sections::Vector{Matrix{Float64}};
		ellipse_points::Int=size(first(sections),1)
	)
		pts = vcat(map(sections) do section
			Point3f.(eachrow(section))
		end...)
		return get_model_contour_mesh(pts; ellipse_points)
	end
	function get_sections(model::AbstractCelegansModel)
		r = LinRange(0, 1, length(model))
		_transverse_splines = transverse_splines(model)
		sections = map(r) do t
			Point3f.(t .|> _transverse_splines)
		end
		return sections
	end
	function get_model_contour_mesh(model::AbstractCelegansModel; kwargs...)
		sections = get_sections(model)
		return get_model_contour_mesh(sections; kwargs...)
	end
	function get_model_contour_mesh(
		sections::Vector{Vector{Point3f}};
		ellipse_points=length(first(sections)),
		kwargs...
	)
		pts = vcat(sections...)
		get_model_contour_mesh(pts; ellipse_points, kwargs...)
	end
	function get_model_contour_mesh(
		pts::Vector{Point3f};
		ellipse_points::Int,
		transform_points = identity
	)
		if transform_points != identity
			pts = transform_points.(pts)
		end
		npts = length(pts)
		#mesh(GeometryBasics.Mesh(pts, faces[1:964]))
		A = mod1.(1:ellipse_points, ellipse_points)
		B = mod1.(2:ellipse_points+1, ellipse_points)
		C = B .+ ellipse_points
		D = A .+ ellipse_points
		faces = QuadFace.(A,B,C,D)
		faces = vec(faces .+
			reshape(QuadFace.(eachrow(
				[ellipse_points ellipse_points ellipse_points ellipse_points] .* (0:length(pts)÷ellipse_points-2)
			)), 1,:)
		)
		M = GeometryBasics.Mesh(pts, faces)
	end
	function get_model_manifold_mesh_components(
		model::AbstractCelegansModel;
		transform_points = identity
	)
		r = LinRange(0, 1, length(model))
		_splines = [transverse_spline(model,1) central_spline(model) transverse_spline(model,17)]
		sections = map(r) do t
			Point3f.(t .|> _splines)
		end
		pts = transform_points.(Iterators.flatten(sections))
		offsets = (0:length(pts)÷3-2) .* (QuadFace(3,3,3,3),)
		faces = QuadFace[(1,2,5,4) (2,3,6,5)] .+ offsets
		return pts, vec(faces)
	end
	function get_model_manifold_mesh(
		model::AbstractCelegansModel;
		transform_points = identity
	)
		pts, faces = get_model_manifold_mesh_components(model; transform_points)
		return GeometryBasics.Mesh(pts, faces)
	end