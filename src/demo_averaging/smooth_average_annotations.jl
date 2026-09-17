include("smooth_polar_dct1.jl")

_endpoint_pin_weight(d::Int, taper_width::Int) =
    d >= taper_width ? 0.0 : (cos(π * d / taper_width) + 1) / 2

# Pins `smoothed` to `raw`'s exact endpoint values by adding a smooth
# (raised-cosine) correction that equals the raw/smoothed mismatch at t=1/t=N
# and decays to 0 by `taper_width` timepoints in. Blending toward `raw` itself
# (as an earlier version of this did) reintroduces raw's full high-frequency
# content everywhere inside the taper window, not just at the two endpoints;
# this correction is a constant offset shaped by a smooth taper, so it carries
# no per-timepoint noise from `raw`.
function _blend_to_raw_endpoints(smoothed, raw, taper_width::Int)
    taper_width <= 0 && return smoothed
    N = length(smoothed)
    offset_start = raw[1] .- smoothed[1]
    offset_end = raw[end] .- smoothed[end]
    return map(eachindex(smoothed)) do t
        w_start = _endpoint_pin_weight(t - 1, taper_width)
        w_end = _endpoint_pin_weight(N - t, taper_width)
        smoothed[t] .+ w_start .* offset_start .+ w_end .* offset_end
    end
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