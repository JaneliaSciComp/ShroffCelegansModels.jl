"""
    MIPAVIO

MIPAV (Medical Image Processing, Analysis, and Visualization) Input / Output
utility module.

MIPAV is a Java program from the National Institutes of Health
"""
module MIPAVIO
    using DataFrames: DataFrame
    using GeometryBasics: Point3
    using CSV: CSV
    using ShroffCelegansModels: Datasets

    export mipav_df_to_points, mipav_df_to_point_dict

    function mipav_df_to_points(df::DataFrame)
        map(df.x_voxels, df.y_voxels, df.z_voxels) do x,y,z
            Point3(x,y,z)
        end
    end

    function mipav_df_to_point_dict(df::DataFrame)
        map(df.name, df.x_voxels, df.y_voxels, df.z_voxels) do name, x,y,z
            name => Point3(x,y,z)
        end |> Dict{eltype(df.name), Point3{eltype(df.x_voxels)}}
    end

    function get_integrated_annotations(ds::Datasets.NormalizedDataset, time_offset=1)
        data_path = joinpath("integrated_annotation","annotations.csv")
        return CSV.read(get_model_csv(ds, data_path, time_offset), DataFrame)
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
            @info "Filepath" filepath
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
end