using Printf
using GeometryBasics
using ProgressMeter
using CSV
using DataFrames

include("get_group_annotation_positions_over_time.jl")

function average_annotations(datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset})
    cache = my_annotation_position_cache
    group_annotation_positions_over_time = get_group_annotation_positions_over_time(datasets, cache)
    #common_annotations = intersect(map(datasets_info) do dataset_info
    #    collect(keys(dataset_info.annotation_dict))
    #end...)
    datasets_info = get_datasets_info(datasets)

    # common_annotations
    annotations = intersect(map(datasets_info) do dataset_info
        collect(keys(dataset_info.annotation_dict))
    end...)

    positions = map(eachindex(first(group_annotation_positions_over_time))) do j
        map(annotations) do name
            mean(map(eachindex(group_annotation_positions_over_time)) do i
                group_annotation_positions_over_time[i][j][name]
            end)
        end
    end
    return (; annotations, positions)
end

# average_annotations_dict = average_annotations(datasets)
function average_annotations(datasets::Dict{String, Vector{ShroffCelegansModels.Datasets.NormalizedDataset}})
    average_annotations_dict = Dict(keys(datasets) .=> map(collect(keys(datasets))) do k    
           average_annotations(datasets[k])
    end)
    return average_annotations_dict
end

# save_average_annotations(average_annotations_dict; filename = "average_annotations.h5")
function save_average_annotations(
    average_annotations_dict::Dict{String, @NamedTuple{annotations::Vector{String}, positions::Vector{Vector{Point{3, Float64}}}}};
    filename = "average_annotations.h5"
)
    h5open(filename, "w") do h5f
        for (k, v) in average_annotations_dict
            h5g = create_group(h5f, k)
            h5g["annotations"] = v.annotations
            for (idx, points) in pairs(v.positions)
                _points = reinterpret(Float64, points)
                _points = reshape(_points, 3, :)
                _points = transpose(_points)
                tp_name = @sprintf("timepoint_%03d", idx)
                h5g[tp_name] = collect(_points)
            end
        end
    end
end

function load_average_annotations(; filename = "average_annotations.h5")
    d = Dict{String, @NamedTuple{annotations::Vector{String}, positions::Vector{Vector{Point{3, Float64}}}}}()
    h5open(filename) do h5f
        for k in keys(h5f)
            h5g = h5f[k]
            annotations = h5g["annotations"][]::Vector{String}
            timepoints = filter(contains("timepoint_"), keys(h5g))
            positions = Vector{Vector{Point3{Float64}}}(undef, length(timepoints))
            for tp in timepoints
                matrix = transpose(h5g[tp][]::Matrix{Float64})
                idx = parse(Int, tp[end-2:end])
                positions[idx] = vec(reinterpret(Point3{Float64}, matrix))
            end
            d[k] = (; annotations, positions)
        end
    end
    return d
end
