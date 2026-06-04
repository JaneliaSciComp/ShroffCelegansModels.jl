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

    # Generic per-timepoint mtime scanner. `path_builder(dataset, time_offset)`
    # returns the file path(s) to stat for that timepoint — `String`, `Vector{String}`,
    # `missing`, or empty. The returned mtime is the *max* across the listed paths
    # (so a dataset's mtime advances when any of its inputs are touched). NaN if
    # nothing stat-able.
    #
    # `path_builder` is the first argument so callers can use Julia do-block syntax:
    #
    #     _mtimes_unix(dataset) do ds, i
    #         get_integrated_annotations_path(ds, i)
    #     end
    function _mtimes_unix(path_builder, dataset::Datasets.NormalizedDataset)::Vector{Float64}
        map(1:length(range(dataset.cell_key))) do i
            try
                paths = path_builder(dataset, i)
                ismissing(paths) && return NaN
                if paths isa AbstractString
                    return stat(paths).mtime
                end
                isempty(paths) && return NaN
                maximum(stat(p).mtime for p in paths)
            catch
                NaN
            end
        end
    end

    # The lattice-related files consumed by `get_avg_models` for a single
    # timepoint: `lattice_final/lattice.csv` plus every
    # `model_crossSections/latticeCrossSection_*.csv`. Returns `missing` for
    # outlier timepoints (no Decon_reg_$tp dir on disk) — same convention as
    # `get_model_csv`. The cross-section dir may not exist on some datasets;
    # in that case we just return the lattice file alone.
    function _lattice_paths(ds::Datasets.NormalizedDataset, time_offset::Int)
        timepoint = range(ds.cell_key)[time_offset]
        timepoint ∈ ds.cell_key.outliers && return missing
        base = joinpath(ds.path, "Decon_reg_$(timepoint)")
        lattice = joinpath(base, "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
        cs_dir = joinpath(base, "model_crossSections")
        cs_paths = if isdir(cs_dir)
            joinpath.(cs_dir,
                      filter(f -> startswith(f, "latticeCrossSection_") && endswith(f, ".csv"),
                             readdir(cs_dir)))
        else
            String[]
        end
        return String[lattice; cs_paths...]
    end

    get_annotation_modified_times_unix(dataset::Datasets.NormalizedDataset) =
        _mtimes_unix(dataset) do d, i
            get_integrated_annotations_path(d, i)
        end

    get_lattice_modified_times_unix(dataset::Datasets.NormalizedDataset) =
        _mtimes_unix(_lattice_paths, dataset)

    function get_annotation_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}}
    )::Dict{String, Vector{Vector{Float64}}}
        Dict(k => get_annotation_modified_times_unix.(v) for (k, v) in datasets)
    end

    function get_lattice_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}}
    )::Dict{String, Vector{Vector{Float64}}}
        Dict(k => get_lattice_modified_times_unix.(v) for (k, v) in datasets)
    end

    # Back-compat alias. Existing callers expecting "the annotation mtimes" keep
    # working unchanged.
    get_modified_times_unix(dataset::Datasets.NormalizedDataset) =
        get_annotation_modified_times_unix(dataset)
    get_modified_times_unix(datasets::Dict{String, Vector{Datasets.NormalizedDataset}}) =
        get_annotation_modified_times_unix(datasets)

    # Write a single kind's mtimes into an HDF5 group named `kind` (e.g.
    # "annotation", "lattice"). Supports both fresh-file ("w") and append ("r+")
    # modes via the `mode` arg so `save_all_modified_times_unix` can stack
    # multiple kinds in one file.
    function save_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}},
        modified_times::Dict{String, Vector{Vector{Float64}}} = get_modified_times_unix(datasets);
        filepath::String,
        kind::String = "annotation",
        mode::String = "w",
    )
        h5open(filepath, mode) do h5f
            kind_group = create_group(h5f, kind)
            for group in keys(datasets)
                h5g = create_group(kind_group, group)
                for (k, v) in pairs(modified_times[group])
                    h5g[string(k)] = v
                    A = attrs(h5g[string(k)])
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

    # Convenience: scan both annotation and lattice mtimes and write them to a
    # single file under top-level groups `annotation/` and `lattice/`. Used by
    # the daily save-modified-times CronJob.
    function save_all_modified_times_unix(
        datasets::Dict{String, Vector{Datasets.NormalizedDataset}};
        filepath::String,
        annotation_mtimes::Dict{String, Vector{Vector{Float64}}} = get_annotation_modified_times_unix(datasets),
        lattice_mtimes::Dict{String, Vector{Vector{Float64}}} = get_lattice_modified_times_unix(datasets),
    )
        save_modified_times_unix(datasets, annotation_mtimes; filepath, kind="annotation", mode="w")
        save_modified_times_unix(datasets, lattice_mtimes;    filepath, kind="lattice",    mode="r+")
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