using ShroffCelegansModels: CelegansModel, Datasets

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

function annotation_positions(smts_nt, annotation_dict, nt; avg_models::Vector{<: CelegansModel} = avg_models, ws=nothing)
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
            end; ws
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
    avg_models::Vector{<: CelegansModel} = avg_models,
    checkpoint_dir::Union{Nothing, AbstractString} = nothing,
    progress_counter::Union{Nothing, Threads.Atomic{Int}} = nothing,
    progress_total::Union{Nothing, Int} = nothing,
)::Vector{Vector{Dict{String, Point3{Float64}}}}
    @assert length(avg_models) == length(normalized_timepoints)

    # `get_datasets_info` runs sequentially before the threaded loop. Time it
    # explicitly so we can see how much of step 6's budget is serial setup
    # (StraightenedModelTimeSeries + get_cell_trajectory_dict per dataset).
    @info "get_datasets_info begin" n_datasets=length(datasets)
    t_info = time()
    datasets_info = get_datasets_info(datasets)
    @info "get_datasets_info done" elapsed_s=round(time() - t_info; digits=2)

    # Global progress numbering: when the caller (the Dict-keyed
    # average_annotations) threads a shared counter + grand total through, each
    # "step6 dataset done" line is numbered against ALL datasets across every
    # group, so overall progress is estimable. Standalone calls fall back to
    # this group's own dataset count.
    grand_total = progress_total === nothing ? length(datasets_info) : progress_total
    done_counter = progress_counter === nothing ? Threads.Atomic{Int}(0) : progress_counter

    # Parallelize across datasets — the inner per-timepoint loop is replaced
    # with a sequential map so we don't oversubscribe threads. Each dataset's
    # work is independent except for the shared `cache` write, which is
    # guarded by a lock.
    out = Vector{Vector{Dict{String, Point3{Float64}}}}(undef, length(datasets_info))
    cache_lock = ReentrantLock()
    prog = ProgressMeter.Progress(length(datasets_info); desc="Avg annotations / dataset...")
    Threads.@threads for ds_idx in eachindex(datasets_info)
        ds_t0 = time()
        dataset_info = datasets_info[ds_idx]
        dataset = dataset_info.dataset
        annotation_dict = dataset_info.annotation_dict
        smts_nt = dataset_info.smts_nt
        cached = lock(cache_lock) do
            get(cache, dataset.path, nothing)
        end
        from_cache = cached !== nothing
        _annotation_positions_over_time = if from_cache
            cached
        else
            local positions = Vector{Vector{Point3{Float64}}}(undef, length(normalized_timepoints))
            # One TPS workspace per dataset (= per @threads task ⇒ thread-local),
            # reused across this dataset's timepoints so the ~tens-of-MB solve
            # buffers aren't re-allocated ~371× here. Lazily sized on first solve.
            local ws = Ref{Any}(nothing)
            for i in eachindex(normalized_timepoints)
                positions[i] = annotation_positions(smts_nt, annotation_dict, normalized_timepoints[i]; avg_models, ws)
            end
            lock(cache_lock) do
                cache[dataset.path] = positions
            end
            positions
        end
        _annotation_positions_over_time::Vector{Vector{Point3{Float64}}}
        out[ds_idx] = try
            map(_annotation_positions_over_time) do points
                Dict{String, Point3d}(keys(annotation_dict) .=> points)
            end
        catch err
            @error "Error creating dict for dataset $(dataset.path)" exception = (err, Base.catch_backtrace())
            map(_annotation_positions_over_time) do points
                Dict{String, Point3d}(keys(annotation_dict) .=> fill(Point3(NaN), length(keys(annotation_dict))))
            end
        end
        # Write per-dataset checkpoint (skip if entry came from a prior checkpoint
        # load — no point rewriting unchanged data). Each thread writes a unique
        # filename so no lock is needed for the file I/O itself.
        checkpoint_bytes = if checkpoint_dir !== nothing && !from_cache
            ShroffCelegansModels.write_dataset_checkpoint(
                checkpoint_dir, dataset.path, _annotation_positions_over_time,
            )
        else
            0
        end
        n_done = Threads.atomic_add!(done_counter, 1) + 1
        @info("step6 dataset done",
            dataset = string(n_done, "/", grand_total),
            group_idx = ds_idx,
            group_total = length(datasets_info),
            path = dataset.path,
            thread = Threads.threadid(),
            n_annotations = length(annotation_dict),
            n_timepoints = length(normalized_timepoints),
            elapsed_s = round(time() - ds_t0; digits=2),
            from_cache = from_cache,
            checkpoint_bytes = checkpoint_bytes,
        )
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)
    return out
end
