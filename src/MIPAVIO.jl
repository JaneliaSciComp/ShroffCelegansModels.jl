"""
    MIPAVIO

MIPAV (Medical Image Processing, Analysis, and Visualization) Input / Output
utility module.

MIPAV is a Java program from the National Institutes of Health
"""
module MIPAVIO
    using DataFrames: DataFrame
    using GeometryBasics: Point3, Point3f, Point3d
    using CSV: CSV
    using ShroffCelegansModels: Datasets
    using Statistics: mean
    using Dates: DateTime, unix2datetime, TimeType
    using HDF5: h5open, create_group, attrs

    export mipav_df_to_points, mipav_df_to_point_dict

    function mipav_df_to_points(df::DataFrame)::Vector{Point3d}
        map(
            df.x_voxels::Vector{Float64},
            df.y_voxels::Vector{Float64},
            df.z_voxels::Vector{Float64}
        ) do x,y,z
            Point3d(x,y,z)
        end
    end

    function mipav_df_to_point_dict(df::DataFrame)
        map(df.name, df.x_voxels, df.y_voxels, df.z_voxels) do name, x,y,z
            name => Point3(x,y,z)
        end |> Dict{eltype(df.name), Point3{eltype(df.x_voxels)}}
    end

    function get_integrated_annotations_path(ds::Datasets.NormalizedDataset, time_offset=1)
        data_path = joinpath("integrated_annotation","annotations.csv")
        return get_model_csv(ds, data_path, time_offset)
    end

    function get_integrated_annotations(ds::Datasets.NormalizedDataset, time_offset=1; validate = true)
        path = get_integrated_annotations_path(ds, time_offset)
        df = CSV.read(path, DataFrame)
        if validate
            names = df.name
            for k in keys(ds.cell_key.mapping)
                if string(k) ∉ names
                    #error("$k is not in $(ds.path)")
                    @warn "$k is not in $(path)"
                end
            end
        end
        return df
    end

    function get_integrated_annotations(::Type{Dict}, args...)
        df = get_integrated_annotations(args...)
        return Dict(row[1] => Point3(row[2], row[3], row[4]) for row in eachrow(Matrix(df)))
    end

    function get_model_csv(ds::Datasets.NormalizedDataset, data_path, time_offset=1)::Union{Missing, String}
        timepoint = range(ds.cell_key)[time_offset]
        if timepoint ∈ ds.cell_key.outliers
            return missing
        else
            filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", data_path)
            @debug "Filepath" filepath
            if isfile(filepath)
                return filepath
            else
                throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
            end
        end
    end

    function get_straightened_lattice(ds::Datasets.NormalizedDataset, time_offset=1)
        data_path = joinpath("straightened_lattice", "straightened_lattice.csv")
        return CSV.read(get_model_csv(ds, data_path, time_offset), DataFrame)
    end

    function get_straightened_lattice_xy_center(ds::Datasets.NormalizedDataset, time_offset=1)
        csv = get_straightened_lattice(ds, time_offset)
        m = Matrix(csv[1:2:end, 2:3] .+ csv[2:2:end, 2:3])./2
        Point3f(mean(eachrow(m))..., 0)
    end

    #=
    function get_straightened_lattice_xy_center(ds::Datasets.NormalizedDataset, time_offset=1)
        csv = get_straightened_lattice(ds, time_offset)
        m = Matrix(csv[1:2:end, 2:3] .+ csv[2:2:end, 2:3])./2
        Point3f(mean(eachrow(m))..., 0)
    end
    =#

    function get_straightened_annotations(ds::Datasets.NormalizedDataset, time_offset=1)::Union{Missing, String}
        timepoint = range(ds.cell_key)[time_offset]
        if timepoint ∈ ds.cell_key.outliers
            return missing
        else
            filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "straightened_annotations", "straightened_annotations.csv")
            if isfile(filepath)
                return filepath
            else
                throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
            end
        end
    end

    function load_straightened_annotations_over_time(dataset::Datasets.NormalizedDataset, offsets = 1:length(range(dataset.cell_key)))
        annotations = map(offsets) do time_offset
            path = get_straightened_annotations(dataset, time_offset)
            if ismissing(path)
                return missing
            end
            annotation_df = CSV.read(path, DataFrame)
            pts = Point3f.(eachrow(Matrix(annotation_df)[:, 2:4]))
            pts .-= get_straightened_lattice_xy_center(dataset, time_offset)
            Dict(annotation_df[:,1] .=> pts)
        end
        return annotations
    end

    function get_modified_times_unix(dataset::Datasets.NormalizedDataset)::Vector{Float64}
        map(1:length(range(dataset.cell_key))) do i
            try
                path = get_integrated_annotations_path(dataset, i)
                ismissing(path) && return NaN
                path_stat = stat(path)
                path_stat.mtime
            catch
                NaN
            end
        end
    end

    function get_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}}
    )::Dict{String,Vector{Vector{Float64}}}
        Dict(k => get_modified_times_unix.(v) for (k,v) in datasets)
    end

    function save_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}},
        modified_times::Dict{String,Vector{Vector{Float64}}} = get_modified_times_unix(datasets);
        filepath::String
    )
        h5open(filepath, "w") do h5f
            for group in keys(datasets)
                h5g = create_group(h5f, group)
                for (k,v) in pairs(modified_times[group])
                    h5g[string(k)] = v
                    A = attrs(h5f[group][string(k)])
                    dataset = datasets[group][k]
                    A["path"] = dataset.path
                    A["cell_key.name"] = dataset.cell_key.name
                    A["cell_key.start"] = dataset.cell_key.start
                    A["cell_key.end"] = dataset.cell_key.stop
                    A["cell_key.outliers"] = dataset.cell_key.outliers
                end
            end
        end
    end

    function get_modified_times(dataset::Datasets.NormalizedDataset)::Vector{Union{Missing,DateTime}}
        map(1:length(range(dataset.cell_key))) do i
            path = get_integrated_annotations_path(dataset, i)
            ismissing(path) && return missing
            path_stat = stat(path)
            unix2datetime(path_stat.mtime)
        end
    end

    function get_last_modified_time(dataset::Datasets.NormalizedDataset)
        mtimes = get_modified_times(dataset)
        if all(ismissing, mtimes)
            return missing
        else
            return maximum(skipmissing(mtimes))
        end
    end

    function print_integrated_annotations_modified_since(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}},
        since::TimeType;
        modified_times = get_modified_times_unix(datasets)
    )
        for group in keys(modified_times)
            for idx in keys(modified_times[group])
                for tp in keys(modified_times[group][idx])
                    ds = datasets[group][idx]
                    u = modified_times[group][idx][tp]
                    if !isnan(u) && unix2datetime(u) > since
                        path = get_integrated_annotations_path(ds, idx)
                        println(group, ", ", idx, ", ", tp, ", ", unix2datetime(u))
                        # println(path)
                    end
                end
            end
        end
    end
end