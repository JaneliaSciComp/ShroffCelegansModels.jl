using ShroffCelegansModels.HDF5
using ShroffCelegansModels.StatsBase: zscore
using ShroffCelegansModels.DataFrames
using ShroffCelegansModels.CSV
using ShroffCelegansModels.GeometryBasics
using Printf

function zscore_analysis()
    df = DataFrame(dataset=String[], embryo=String[], annotation=String[], zscore=Float64[])
    h5open("embryos_371_2026_04_10.h5") do h5f
        for k in keys(h5f)
            for e in keys(h5f[k])
                M = h5f[k][e][3, :, :]
                zscores = zscore(sqrt.(sum(diff(M; dims=1) .^ 2; dims=1)))
                pairs = Dict(attrs(h5f[k][e])["annotation_names"] .=> zscores')
                deviation = filter(k -> pairs[k] > 2, keys(pairs))
                for d in deviation
                    println(k, ", ", e, ", ", d, ", ", pairs[d])
                    push!(df, (k, e, d, pairs[d]))
                end
            end
        end
    end
    return df
end

function get_annotations_cache_keys(datasets::Dict{String, Dict{String, ShroffCelegansModels.Dataset}})
    flatened_datasets = [dataset for group in values(datasets) for dataset in values(group)]
    dataset_keys = ShroffCelegansModels.annotations_cache_key.(flattened_datasets)
end
function raw_zscore_analysis(
    datasets::Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}};
    threshold=1,
    time_threshold=1
)
    d = map(collect(keys(datasets))) do group
        group => map(collect(keys(datasets[group]))) do embryo
            dataset = datasets[group][embryo]
            dict = raw_annotation_dict(dataset)
            map(collect(keys(dict))) do annotation
                timeseries = dict[annotation]
                annotation => sqrt(sum(diff((x->x[3]).(skipmissing(timeseries))).^2))
            end |> Dict
        end
    end |> Dict
    flattened_dict = Dict(
        (group,embyro,annotation) => value
        for (group, embyros) in d
            for (embyro,annotations) in pairs(embyros) 
                for (annotation,value) in annotations
    )
    zscore_dict = Dict(keys(flattened_dict) .=> zscore(values(flattened_dict) |> collect))
    filtered_zscore_dict = filter(zscore_dict) do (k, zscore)
        zscore > threshold
    end
    map(keys(filtered_zscore_dict) |> collect) do k
        (group, embryo, annotation) = k
        dataset = datasets[group][embryo]
        dict = raw_annotation_dict(dataset)
        timeseries = dict[annotation]
        timepoint = find_max_diff_idx(timeseries) + dataset.cell_key.start - 1
        # println(group, ", ", embryo, ", ", annotation, ", ", zscore, ", ", timepoint)
        # get_fix_url(group, embryo, annotation, timepoint)
        (;group, embryo, annotation, timepoint, zscore = zscore_dict[k], link = get_fix_url(group, embryo, annotation, timepoint))
    end |> DataFrame
end
function raw_annotation_dict(dataset; use_myuntwist = true)
    straighted_annotations_over_time = ShroffCelegansModels.load_straightened_annotations_over_time(dataset; use_myuntwist)
    mapping = dataset.cell_key.mapping
    timeseries = map(keys(mapping) |> collect) do annotation_name
        Vector{Union{Missing, Point3d}}(map(straighted_annotations_over_time) do d
            if ismissing(d)
                return missing
            end
            get(d, string(annotation_name), missing)
        end)
    end::Vector{Vector{Union{Missing, Point3d}}}

    return Dict(values(mapping) .=> timeseries)
end
function find_max_diff_idx(timeseries::Vector{Union{Missing, Point3d}})
    idx = [i for (i, point) in enumerate(timeseries) if !ismissing(point)]
    delta = diff((x->x[3]).(skipmissing(timeseries)))
    abs_delta = abs.(delta)
    i_max = argmax(abs_delta)
    if i_max == 1 || i_max == length(idx)
        return idx[i_max]
    end
    if delta[i_max] > 0
        if delta[i_max-1] < -delta[i_max]/2
            return idx[i_max]
        else
            return idx[i_max+1]
        end
    else
        if delta[i_max-1] > -delta[i_max]/2
            return idx[i_max]
        else
            return idx[i_max+1]
        end
    end
end
function find_max_diff_idx(
    dataset::ShroffCelegansModels.NormalizedDataset,
    annotation::String,
    dict::Dict{String, Vector{Union{Missing, Point3d}}} = raw_annotation_dict(dataset)
)
    find_max_diff_idx(dict[annotation]) + dataset.cell_key.start - 1
end
function zscore_dict(dataset; use_myuntwist = true)
    dict = raw_annotation_dict(dataset; use_myuntwist)
    zscores = zscore(map(values(dict)) do timeseries
        sqrt(sum(diff((x->x[3]).(skipmissing(timeseries))).^2))
    end)
    indices = map(values(dict)) do timeseries
        find_max_diff_idx(timeseries) + dataset.cell_key.start - 1
    end

    return Dict(keys(dict) .=> zip(zscores, indices))
end
function outlier_dict(dataset; use_myuntwist = true)
    dict = raw_annotation_dict(dataset; use_myuntwist)
    indices = map(values(dict)) do timeseries
        find_max_diff_idx(timeseries) + dataset.cell_key.start - 1
    end
    return Dict(keys(dict) .=> indices)
end
function get_fix_url(group, group_idx, annotation, timepoint)
    return @sprintf(
        "https://%s/fix_annotation_ap_axis/%s/%d?annotation=%s&timepoint=%d",
        get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"),
        group,
        group_idx,
        annotation,
        timepoint
    )
end
function dataset_to_links(datasets, group, group_idx, threshold; use_myuntwist = true)
    dataset = datasets[group][group_idx]
    dict = zscore_dict(dataset; use_myuntwist)
    filter!(dict) do (annotation, (zscore, _))
        zscore > threshold
    end
    links = map(collect(keys(dict))) do annotation
        (zscore, timepoint) = dict[annotation]
        get_fix_url(group, group_idx, annotation, timepoint)
    end
    return links
end
function list_outliers(;threshold=1)
    for k1 in keys(datasets)
        for k2 in keys(datasets[k1])
            dataset_to_links(datasets, k1, k2, threshold) .|> println;
        end
    end
end
