using Printf
using GeometryBasics
using ProgressMeter
using CSV
using DataFrames
using HDF5

# include("get_group_annotation_positions_over_time.jl")
using ShroffCelegansModels: CelegansModel, get_datasets_info, get_group_annotation_positions_over_time, annotation_positions

function average_annotations(
    datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset};
    cache::Dict{String, Vector{Vector{Point3{Float64}}}} = my_annotation_position_cache,
    timepoints::Union{AbstractVector{Float64}, Integer} = LinRange(0,1,201),
    avg_models::Vector{<: CelegansModel} = avg_models,
    use_cell_key_annotations_only = true,
    checkpoint_dir::Union{Nothing, AbstractString} = nothing,
)
    if isa(timepoints, Integer)
        N_timepoints = timepoints
        timepoints = LinRange(0, 1, N_timepoints)
    end
    group_annotation_positions_over_time = get_group_annotation_positions_over_time(
        datasets, cache, timepoints;
        avg_models = avg_models,
        checkpoint_dir = checkpoint_dir,
    )
    group_annotation_positions_over_time::Vector{Vector{Dict{String, Point3{Float64}}}}
    #common_annotations = intersect(map(datasets_info) do dataset_info
    #    collect(keys(dataset_info.annotation_dict))
    #end...)
    datasets_info = get_datasets_info(datasets)

    # common_annotations
    annotations = intersect(map(datasets_info) do dataset_info
        dataset_annotations = collect(keys(dataset_info.annotation_dict))
        if use_cell_key_annotations_only
            dataset_annotations = filter(name -> name ∈ values(dataset_info.dataset.cell_key.mapping), dataset_annotations)
        end
        dataset_annotations
    end...)::Vector{String}

    positions = map(eachindex(first(group_annotation_positions_over_time))) do j
        map(annotations) do name
            mean(map(eachindex(group_annotation_positions_over_time)) do i
                group_annotation_positions_over_time[i][j][name]
            end)
        end::Vector{Point3{Float64}}
    end::Vector{Vector{Point3{Float64}}}
    return (; annotations, positions)
end

# average_annotations_dict = average_annotations(datasets)
function average_annotations(
    datasets::Dict{String, Vector{ShroffCelegansModels.Datasets.NormalizedDataset}};
    cache::Dict{String, Vector{Vector{Point3{Float64}}}} = my_annotation_position_cache,
    timepoints::Union{AbstractVector{Float64}, Integer} = LinRange(0,1,201),
    avg_models::Vector{<: CelegansModel} = avg_models,
    use_cell_key_annotations_only = true,
    checkpoint_dir::Union{Nothing, AbstractString} = nothing,
)
    average_annotations_dict = Dict(keys(datasets) .=> map(collect(keys(datasets))) do k
           average_annotations(datasets[k]; cache, timepoints, avg_models, use_cell_key_annotations_only, checkpoint_dir)
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

"""
    load_latest_average_annotations(; default_filename, prefix="edited_smoothed_average_annotations_", dir=ENV["RECOMPUTE_OUTPUT_DIR"] or "/data/annotations/recompute")

Load the most recently produced averaged-annotations HDF5 from the recompute
output directory if one exists; otherwise fall back to `default_filename`.
Used by the meshscatter web apps so a fresh recompute is picked up on the
next service restart without code changes.

The directory is scanned for files matching `prefix*.h5` and the newest by
mtime is chosen.
"""
function load_latest_average_annotations(;
    default_filename::AbstractString,
    prefix::AbstractString = "edited_smoothed_average_annotations_",
    dir::AbstractString = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute"),
)
    chosen = default_filename
    if isdir(dir)
        candidates = String[
            joinpath(dir, f) for f in readdir(dir)
            if startswith(f, prefix) && endswith(f, ".h5") && !occursin(".tmp.", f)
        ]
        if !isempty(candidates)
            chosen = argmax(mtime, candidates)
        end
    end
    @info "Loading averaged annotations" chosen default_filename dir
    return load_average_annotations(; filename = chosen)
end
