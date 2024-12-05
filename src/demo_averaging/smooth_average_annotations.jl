include("smooth_polar_dct1.jl")

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
    smooth_factor_z = 0.04
)
    smoothed = DICT_TYPE()
    for (k,v) in average_annotations_dict
        smoothed_positions = map(v.positions |> stack |> eachrow) do positions_over_time
            smooth_polar_dct1(
                positions_over_time,
                1/smooth_factor_r,
                1/smooth_factor_θ,
                1/smooth_factor_z
            )
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