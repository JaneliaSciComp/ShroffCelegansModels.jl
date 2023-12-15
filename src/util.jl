const get_circle_points = let N=32, θ = (0:N-1)//N*2, cθ = cospi.(θ), sθ = sinpi.(θ)
	function _get_circle_points(right_vector, normal_vector, center_point)
		# This could be precomputed once
		#=
		θ = LinRange(0, 2π, 33)[1:end-1]
		cθ = cos.(θ)
		sθ = sin.(θ)
		=#
		map(cθ, sθ) do c, s
			# This is done to match mipav
			right_vector*c + normal_vector*s + center_point
		end
	end
end