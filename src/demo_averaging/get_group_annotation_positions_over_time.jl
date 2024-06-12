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
                title_twisted[] = "Twisted; idx = $nt"
                mts(nt)
            end
        end


        annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)
        _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
        annotation_dict = Dict(_annotation_text .=> values(annotation_dict))
        return (; dataset, smts_nt, mts_nt, annotation_dict)
    end
end
function get_group_annotation_positions_over_time(datasets, cache)
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
                title_twisted[] = "Twisted; idx = $nt"
                mts(nt)
            end
        end


        annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)
        _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
        annotation_dict = Dict(_annotation_text .=> values(annotation_dict))
        return (; dataset, smts_nt, mts_nt, annotation_dict)
    end
    group_annotation_positions_over_time = map(datasets_info) do dataset_info
        dataset = dataset_info.dataset
        annotation_dict = dataset_info.annotation_dict
        smts_nt = dataset_info.smts_nt
        if haskey(cache, dataset.path)
            _annotation_positions_over_time = cache[dataset.path]
        else
            #_annotation_positions_over_time = annotation_positions.((smts_nt,), (annotation_dict,), r)
            _annotation_positions_over_time = Vector{Vector{Point3{Float64}}}(undef, length(r))
            @showprogress Threads.@threads for i in eachindex(r)
                nt = r[i]
                _annotation_positions_over_time[i] = annotation_positions(smts_nt, annotation_dict, nt)
            end
            cache[dataset.path] = _annotation_positions_over_time
        end
        _annotation_positions_over_time # Vector{Vector{Point3{Float64}}}
        map(_annotation_positions_over_time) do positions
            Dict(keys(annotation_dict) .=> positions)
        end
    end # Vector{Vector{Dict{String, Point3{Float64}}}} # dataset, normalized time, name => position
end