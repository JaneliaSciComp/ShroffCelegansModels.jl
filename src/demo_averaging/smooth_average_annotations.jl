include("smooth_polar_dct1.jl")

function _endpoint_taper_weights(N::Int, taper_width::Int)
    w = zeros(Float64, N)
    taper_width <= 0 && return w
    for t in 1:N
        d = min(t - 1, N - t)   # distance from the nearer end, 0 at the very ends
        w[t] = d >= taper_width ? 0.0 : (cos(π * d / taper_width) + 1) / 2
    end
    return w
end

function _blend_to_raw_endpoints(smoothed, raw, taper_width::Int)
    taper_width <= 0 && return smoothed
    w = _endpoint_taper_weights(length(smoothed), taper_width)
    return [w[t] .* raw[t] .+ (1 - w[t]) .* smoothed[t] for t in eachindex(smoothed)]
end

let DICT_TYPE = Dict{
    String,
    @NamedTuple{
        annotations::Vector{String},
        positions::Vector{
            Vector{
                Point{3, Float64}
            }
        }
    }
}

global smooth_average_annotations
function smooth_average_annotations(
    average_annotations_dict::DICT_TYPE;
    smooth_factor_r = 0.05,
    smooth_factor_θ = 0.07,
    smooth_factor_z = 0.04,
    edge_pad::Int = 0,
    taper_width::Int = 0,
)
    smoothed = DICT_TYPE()
    for (k,v) in average_annotations_dict
        smoothed_positions = map(v.positions |> stack |> eachrow) do positions_over_time
            s = smooth_polar_dct1(
                positions_over_time,
                1/smooth_factor_r,
                1/smooth_factor_θ,
                1/smooth_factor_z;
                edge_pad,
            )
            _blend_to_raw_endpoints(s, positions_over_time, taper_width)
        end
        smoothed[k] = (;
            annotations = v.annotations,
            positions = smoothed_positions |> stack |> eachrow |> collect .|> collect
        )
    end
    return smoothed::DICT_TYPE
end

#=
```
smoothed_average_annotations = smooth_average_annotations(average_annotations_dict)
save_average_annotations(smoothed_average_annotations; filename = "smoothed_average_annotations.h5")
```
=#

end # let