function check_annotation_within_contour(io::IO, dataset; use_myuntwist = true)
    mts = ShroffCelegansModels.ModelTimeSeries(dataset)
    _length = length(range(dataset.cell_key))
    expansion_factor = 1.2
    output_lines = String[]
    println(io, join(["Timepoint", "Key", "Name", "Ratio", "Delta", "Filepath"], ", "))
    for t in 1:_length
        model = mts(t)
        annotations_path = ShroffCelegansModels.MIPAVIO.get_integrated_annotations_path(dataset, t)
        annotation_dict = ShroffCelegansModels.twisted_annotations(dataset, t)
        if ismissing(annotation_dict)
            continue
        end
        pts = collect(values(annotation_dict))
        z, central_pts, dist = ShroffCelegansModels.nearest_central_pt(model, pts, expansion_factor)
        z2, central_pts2, dist2 = ShroffCelegansModels.nearest_central_pt(model, pts, 2.0)
        map(zip(keys(annotation_dict), z2, dist, dist2)) do (key, z2, dist, dist2)
            if !isfinite(dist)
                max_radius = ShroffCelegansModels.max_radius_function(model)(z2)
                ratio = dist2 / max_radius
                dist_delta = dist2 - max_radius
                name = get(dataset.cell_key.mapping, Symbol(key), key)
                #@info "Not matched to central point" t key name ratio annotations_path
                println(io, join(string.([t, key, name, ratio, dist_delta, annotations_path]), ","))
            end
        end
    end
    return nothing
end
function check_annotation_within_contour(filename::AbstractString, dataset; use_myuntwist=true)
    open(filename, "w") do io
        check_annotation_within_contour(io, dataset; use_myuntwist)
    end
end
function check_annotation_within_contour(v::Vector{ShroffCelegansModels.Datasets.NormalizedDataset}; use_myuntwist=true)
    foreach(v) do ds
        try
            p = splitpath(ds.path)
            filename = "out_of_contour_" * join(p[3:end], "_") * ".csv"
            check_annotation_within_contour(filename, ds; use_myuntwist)
        catch err
            @error ds err
        end
    end
end