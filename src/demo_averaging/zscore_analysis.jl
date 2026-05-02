using HDF5
using StatsBase: zscore
using DataFrames
using CSV
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
    for group in keys(datasets)
        for embryo in keys(datasets[group])
            dataset = datasets[group][embryo]
            dict = raw_annotation_dict(dataset)
            zscores = zscore(map(values(dict)) do timeseries
                sqrt(sum(diff((x->x[3]).(skipmissing(timeseries))).^2))
            end)
            annotation_names = keys(dict)
            for (annotation, zscore) in zip(annotation_names, zscores)
                if zscore > threshold
                    #println(group, ", ", embryo, ", ", annotation, ", ", zscore)
                    dataset_to_links(datasets, group, embryo, time_threshold) .|> println
                end
            end
        end
    end
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
function zscore_dict(dataset; use_myuntwist = true)
    dict = raw_annotation_dict(dataset; use_myuntwist)
    zscores = zscore(map(values(dict)) do timeseries
        sqrt(sum(diff((x->x[3]).(skipmissing(timeseries))).^2))
    end)
    indices = map(values(dict)) do timeseries
        idx = [i for (i, point) in enumerate(timeseries) if !ismissing(point)]
        delta = diff((x->x[3]).(skipmissing(timeseries)))
        abs_delta = abs.(delta)
        i_max = argmax(abs_delta)
        #return delta[i_max-1:i_max+1]
        if i_max == 1 || i_max == length(idx)
            return idx[i_max] + dataset.cell_key.start - 1
        end
        if delta[i_max] > 0
            if delta[i_max-1] < -delta[i_max]/2
                return idx[i_max] + dataset.cell_key.start - 1 
            else
                return idx[i_max+1] + dataset.cell_key.start - 1
            end
        else
            if delta[i_max-1] > -delta[i_max]/2
                return idx[i_max] + dataset.cell_key.start - 1
            else
                return idx[i_max+1] + dataset.cell_key.start - 1
            end
        end
    end

    return Dict(keys(dict) .=> zip(zscores, indices))
end
function outlier_dict(dataset; use_myuntwist = true)
    dict = raw_annotation_dict(dataset; use_myuntwist)
    indices = map(values(dict)) do timeseries
        idx = [i for (i, point) in enumerate(timeseries) if !ismissing(point)]
        abs_delta = abs.(diff((x->x[3]).(skipmissing(timeseries))))
        i_max = argmax(abs_delta)
        return idx[i_max]
    end
    return Dict(keys(dict) .=> indices)
end
const base_url = "https://shroff-data.int.janelia.org/fix_annotation_ap_axis/%s/%d?annotation=%s&timepoint=%d"
function get_fix_url(group, group_idx, annotation, timepoint)
    return @sprintf(
        "https://shroff-data.int.janelia.org/fix_annotation_ap_axis/%s/%d?annotation=%s&timepoint=%d",
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