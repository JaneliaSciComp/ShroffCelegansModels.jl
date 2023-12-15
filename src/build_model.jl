# ╔═╡ cc0554ba-27b7-4669-8053-be473cf96bfa
function build_celegans_model(
	path::String;
	update_cross_sections::Bool = true, 
)
	df_lattice = CSV.read(path, DataFrame)
	names = df_lattice[:,1]
	if update_cross_sections
		cross_section_path = joinpath(dirname(path), "..", "model_crossSections")
		if isdir(cross_section_path)
			cross_section_files = readdir(cross_section_path)
			n_cross_sections = size(df_lattice,1) ÷ 2
			cross_sections = map(0:n_cross_sections-1) do i
				filename = "latticeCrossSection_$i.csv"
				if filename in cross_section_files
					CSV.read(joinpath(cross_section_path, filename), DataFrame; header=2)
				else
					nothing
				end
			end
			return build_celegans_model(df_lattice; cross_sections, names)
		end
	end
	return build_celegans_model(df_lattice; names)
end

function build_celegans_model(df_lattice::DataFrame; cross_sections::Vector{Union{Nothing, DataFrame}} = Union{Nothing, DataFrame}[], names = String3[])
	seam_cell_data = df_lattice[:,2:4]
	left_seam_cells = seam_cell_data[1:2:end, :]
	right_seam_cells = seam_cell_data[2:2:end, :]
	return build_celegans_model(left_seam_cells, right_seam_cells, cross_sections, names)
end

# ╔═╡ 322eb672-bde2-4857-a275-32639f39408a
function build_celegans_model(left_seam_cells, right_seam_cells, cross_sections::Vector=Union{Nothing,DataFrame}[], names::Vector{<:AbstractString} = String3[])
	P = Point3
	# left_seam_cells = Float32.(left_seam_cells)
	# right_seam_cells = Float32.(right_seam_cells)
	centers = map(eachrow((left_seam_cells .+ right_seam_cells)./2)) do row
		P(row...)
	end
	deltaLengths = norm.(diff(centers))
	afTime = [+0; cumsum(deltaLengths)]
	afTime ./= afTime[end]
	# afTime = Float32.(afTime)

	right_spline = 
		interpolate_natural_cubic_spline(afTime, right_seam_cells)
	left_spline =
		interpolate_natural_cubic_spline(afTime, left_seam_cells)
	center_spline =
		interpolate_natural_cubic_spline(afTime, centers)
	central_spline_d1 = Derivative(1) * center_spline

	# Calculate vectors at each afTime
	forward_vector_afTime = central_spline_d1.(afTime)
	right_vector_afTime = (right_spline.(afTime) - left_spline.(afTime))./2
	radii = norm.(right_vector_afTime)
	normal_vector_afTime = normalize.(cross.(
		forward_vector_afTime,
		right_vector_afTime
	)) .* radii
	circle_points = get_circle_points.(right_vector_afTime, normal_vector_afTime, centers)
	if !isempty(cross_sections)
		CP = eltype(eltype(circle_points))
		for i in eachindex(circle_points)
			if !isnothing(cross_sections[i])
				circle_points[i] = map(eachrow(cross_sections[i])) do row
					CP(row...) + CP(centers[i,:]...)
				end
			end
		end
	end
	n_circle_points = length(first(circle_points))
	contourSplines = map(1:n_circle_points) do c
		contourPoints = map(circle_points) do pts
			pts[c]
		end
		contourSpline = interpolate_natural_cubic_spline(afTime, contourPoints)
	end
	model = CelegansModel(contourSplines, center_spline, names)
	#return contourSplines
	return model
end

build_celegans_model(::Missing) = missing