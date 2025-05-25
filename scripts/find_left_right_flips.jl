function find_left_right_flips(annotation_position_cache = my_annotation_position_cache)
    for (k,v) in my_annotation_position_cache
        first_tp = first(v)::Vector{Point3{Float64}}
        for tp in v
            x = first.(tp)
        end
    end
end

function find_left_right_flips(x::Vector{Float64})
    sx = sign.(x)
    dx = diff(x)
    adx = abs.(dx)
end