function load_straightened_annotations_over_time(dataset::ShroffCelegansModels.Datasets.NormalizedDataset, offsets::UnitRange = 1:length(range(dataset.cell_key)); use_myuntwist::Bool = false)
    key = (dataset.path, offsets, use_myuntwist)
    if haskey(annotations_cache, key)
        return annotations_cache[key]
    end
    if use_myuntwist
        annotations = map(offsets) do time_offset
            ShroffCelegansModels.untwist_annotations(dataset, time_offset)
        end
        annotations_cache[key] = annotations
        return annotations
    else
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
        annotations_cache[key] = annotations
        return annotations
    end
end