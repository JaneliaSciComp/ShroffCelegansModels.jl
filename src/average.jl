using Statistics
using StatsBase
using Missings

const psuedo_seam_cells = ("a0L", "H0L", "H1L", "H2L", "V1L", "V2L", "V3L", "V4L", "V5L", "V6L", "TL")

function average(models::Vector{<: AbstractCelegansModel}, weights::AbstractWeights = uweights(length(models)); n_upsample::Int = 0)
    _names = map(models) do model
        #parent(model).names[1:2:end]
        ShroffCelegansModels.Types.names(model)[1:2:end]
    end
    
    common_names = intersect(_names..., psuedo_seam_cells)
    unique_names = unique(vcat(_names...))
    distinct_names = setdiff(unique_names, common_names)
    #filter!(!startswith("a"), distinct_names)
    # @info "Names" distinct_names
    # common_names = psuedo_seam_cells
    # common_names = distinct_names
    # common_names = intersect(_names...)

    original_knot_pairs = map(models) do model
        twisted_model = parent(model)
        original_knots = knots(central_spline(twisted_model))[4:end-3]
        n = unique(twisted_model.names[1:2:end])
        common_pairs = Iterators.filter(t->(t[2] ∈ common_names),zip(original_knots,n))
        common_names = last.(common_pairs)
        #=
        if length(common_names) > length(psuedo_seam_cells)
            println.(common_names)
        end
        =#
        last.(common_pairs) .=> first.(common_pairs)
    end

    original_knots = map(x->last.(x), original_knot_pairs)

    for i in 1:n_upsample
        original_knots = map(upsample, original_knots)
    end

    #display(original_knots)

    z_positions = map(models, original_knots) do model, original_knots
        # twisted_model = parent(model)
        # original_knots = knots(central_spline(twisted_model))[4:end-3]
        central_spline(model).(original_knots) .|> x->x[3]
    end

    # return z_positions

    # avg_z_position = (z_positions[1] + z_positions[2])./2
    # avg_z_position = sum(z_positions) ./ length(z_positions)
    avg_z_position = mean(z_positions, weights)
    # std_z_position = std(z_positions, weights)

    t_positions = map(models, original_knots) do model, original_knots
        # twisted_model = parent(model)
        # original_knots = knots(central_spline(twisted_model))[4:end-3]
        _transverse_splines = transverse_splines(model)

        map(_transverse_splines) do s
            s.(original_knots)
        end
    end

    t_positions = mean(stack.(t_positions), weights)

    # return t_positions
    new_knots = avg_z_position ./ avg_z_position[end]

    normalized_knots = map(original_knots) do original_knots
        original_knots ./ original_knots[end]
    end

    avg_knots = mean(normalized_knots, weights)
    
    new_knots = avg_knots

    avg_transverse_splines = map(eachcol(t_positions)) do avg_transverse_points
        interpolate_natural_cubic_spline(new_knots, collect(avg_transverse_points))
    end

    avg_central_points = map(avg_z_position) do azp
        Point3(0,0,azp)
    end
    avg_central_spline = interpolate_natural_cubic_spline(new_knots, avg_central_points)

    # return avg_z_position ./ avg_z_position[end]

    # TODO: Change add R names back
    model_names = vec(permutedims([[common_names...] [common_names...]], (2,1)))
    return CelegansModel(avg_transverse_splines, avg_central_spline, model_names)

    # return avg_z_position, t_positions, avg_transverse_splines
    return avg_transverse_splines, avg_central_spline, common_names
end

function upsample(v::Vector{T}) where T
    N = length(v)
    upsampled = T[]
    sizehint!(upsampled, 2N-1)
    for (a,b) in zip(@view(v[1:end-1]), @view(v[2:end]))
        push!(upsampled, a)
        push!(upsampled, (a+b)/2)
    end
    push!(upsampled, v[end])
    return upsampled
end