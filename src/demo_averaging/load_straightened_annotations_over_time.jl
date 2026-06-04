function load_straightened_annotations_over_time(
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset,
    offsets::UnitRange{Int} = 1:length(range(dataset.cell_key));
    use_myuntwist::Bool = false
)::Vector{Union{Missing, Dict{String, Point3d}}}
    key = (dataset.path, offsets, use_myuntwist)
    if haskey(annotations_cache, key)
        return annotations_cache[key].annotations
    end
    if use_myuntwist
        annotations = map(offsets) do time_offset
            ShroffCelegansModels.untwist_annotations(dataset, time_offset)
        end
        annotations_cache[key] = AnnotationsCacheValue(annotations, _dataset_mtime(dataset, offsets))
        return annotations
    else
        annotations = map(offsets) do time_offset
            path = get_straightened_annotations(dataset, time_offset)
            if ismissing(path)
                return missing
            end
            annotation_df = CSV.read(path, DataFrame)
            pts = Point3d.(eachrow(Matrix(annotation_df)[:, 2:4]))
            pts .-= get_straightened_lattice_xy_center(dataset, time_offset)
            Dict{String, Point3d}(annotation_df[:,1] .=> pts)
        end
        annotations_cache[key] = AnnotationsCacheValue(annotations, _dataset_mtime(dataset, offsets))
        return annotations
    end
end

# Dataset-level mtime = max unix mtime across the requested timepoint offsets
# from the integrated_annotation CSVs. NaN if every timepoint is missing.
function _dataset_mtime(
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset,
    offsets::AbstractUnitRange{Int}
)::Float64
    all_mtimes = ShroffCelegansModels.MIPAVIO.get_modified_times_unix(dataset)
    selected = @view all_mtimes[offsets]
    finite = filter(!isnan, selected)
    isempty(finite) ? NaN : maximum(finite)
end

function annotations_cache_key(
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset,
    offsets = 1:length(range(dataset.cell_key)),
    use_myuntwist = true
)
    return (dataset.path, offsets, use_myuntwist)
end