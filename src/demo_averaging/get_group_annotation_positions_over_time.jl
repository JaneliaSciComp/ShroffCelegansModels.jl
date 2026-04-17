include("transform_annotations.jl")

function get_datasets_info(datasets)
    use_myuntwist = true
    datasets_info = map(datasets) do dataset
        smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
        smts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                smts(nt, 2)
            end
        end

        mts = smts.modelTimeSeries
        mts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                nt = round(Int, nt)
                mts(nt)
            end
        end


        annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)
        _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
        annotation_dict = Dict(_annotation_text .=> values(annotation_dict))
        return @NamedTuple{
            dataset::ShroffCelegansModels.Datasets.NormalizedDataset,
            smts_nt::Function,
            mts_nt::Function,
            annotation_dict::Dict{String, Union{Missing, BSplineKit.SplineExtrapolation{<: Spline{Point3d}}}}
        }((dataset, smts_nt, mts_nt, annotation_dict))
    end

    return datasets_info
end

function annotation_positions(smts_nt, annotation_dict, nt; avg_models::Vector{<: CelegansModel} = avg_models)
    N_timepoints = length(avg_models) - 1
    _smodel = smts_nt(nt)
    idx = round(Int, nt * N_timepoints + 1)
    _model = avg_models[idx]
    @debug "annotation positions" nt
    positions = missing
    try
        positions = swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                if ismissing(ann)
                    return Point3(NaN)
                else
                    ann(nt)
                end
            end
        ))
    catch err
        @error "A problem occured at $nt with avg_model[$idx]" exception = (err, Base.catch_backtrace())
        positions = [Point3(NaN) for i in eachindex(annotation_dict)]
    end
    return positions
end

function get_group_annotation_positions_over_time(
    datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset},
    cache::Dict{String, Vector{Vector{Point3{Float64}}}}, # my_annotation_position_cache
    normalized_timepoints::AbstractVector{Float64} = LinRange(0,1,201);
    avg_models::Vector{<: CelegansModel} = avg_models
)::Vector{Vector{Dict{String, Point3{Float64}}}} 
    @assert length(avg_models) == length(normalized_timepoints)
    datasets_info = get_datasets_info(datasets)
    group_annotation_positions_over_time = map(datasets_info) do dataset_info
        dataset = dataset_info.dataset
        annotation_dict = dataset_info.annotation_dict
        smts_nt = dataset_info.smts_nt
        if haskey(cache, dataset.path)
            _annotation_positions_over_time = cache[dataset.path]
        else
            #_annotation_positions_over_time = annotation_positions.((smts_nt,), (annotation_dict,), r)
            _annotation_positions_over_time = Vector{Vector{Point3{Float64}}}(undef, length(normalized_timepoints))
            @showprogress Threads.@threads for i in eachindex(normalized_timepoints)
                nt = normalized_timepoints[i]
                _annotation_positions_over_time[i] = annotation_positions(smts_nt, annotation_dict, nt; avg_models)
            end
            cache[dataset.path] = _annotation_positions_over_time
        end
        _annotation_positions_over_time::Vector{Vector{Point3{Float64}}}
        try
            map(_annotation_positions_over_time) do positions
                Dict{String, Point3d}(keys(annotation_dict) .=> positions)
            end
        catch err
            @error "Error creating dict for dataset $(dataset.path)" exception = (err, Base.catch_backtrace())
            map(_annotation_positions_over_time) do positions
                Dict{String, Point3d}(keys(annotation_dict) .=> fill(Point3(NaN), length(keys(annotation_dict)) ))
            end
        end
    end::Vector{Vector{Dict{String, Point3{Float64}}}} # dataset, normalized time, name => position
    return group_annotation_positions_over_time
end
