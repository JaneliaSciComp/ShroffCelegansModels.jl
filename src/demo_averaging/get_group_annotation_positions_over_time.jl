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

# Compute (or fetch from `cache`) one dataset's positions-over-time. This is the
# per-dataset unit of work shared by the flat `prime_dataset_positions!` pass and
# by `get_group_annotation_positions_over_time`'s loop, so the warp / cache-store /
# checkpoint / log+count logic lives in exactly one place. Returns the raw
# positions vector (the dict-wrapping for the averaging stays in the group caller).
#
# `count_done=false` is used by the post-prime group pass: on a cache hit it logs at
# @debug and does NOT increment the counter, so a primed dataset isn't counted twice
# (prime owns the count). With `count_done=true` (the prime pass and standalone
# callers) every dataset is logged + counted exactly as the original loop did.
function compute_dataset_positions!(
    dataset_info,
    cache::Dict{String, Vector{Vector{Point3{Float64}}}},
    cache_lock::ReentrantLock,
    normalized_timepoints::AbstractVector{Float64},
    avg_models::Vector{<: CelegansModel},
    checkpoint_dir::Union{Nothing, AbstractString};
    done_counter::Threads.Atomic{Int},
    grand_total::Int,
    log_idx::Int = 0,
    log_total::Int = grand_total,
    count_done::Bool = true,
)::Vector{Vector{Point3{Float64}}}
    ds_t0 = time()
    dataset = dataset_info.dataset
    annotation_dict = dataset_info.annotation_dict
    smts_nt = dataset_info.smts_nt

    cached = lock(cache_lock) do
        get(cache, dataset.path, nothing)
    end
    from_cache = cached !== nothing
    positions = if from_cache
        cached
    else
        local pos = Vector{Vector{Point3{Float64}}}(undef, length(normalized_timepoints))
        # One TPS workspace per dataset (= per @threads task ⇒ thread-local),
        # reused across this dataset's timepoints so the ~tens-of-MB solve buffers
        # aren't re-allocated ~371× here. Lazily sized on first solve.
        local ws = Ref{Any}(nothing)
        for i in eachindex(normalized_timepoints)
            pos[i] = annotation_positions(smts_nt, annotation_dict, normalized_timepoints[i]; avg_models, ws)
        end
        lock(cache_lock) do
            cache[dataset.path] = pos
        end
        pos
    end
    positions::Vector{Vector{Point3{Float64}}}

    # Write per-dataset checkpoint (skip cache hits — no point rewriting unchanged
    # data). Each thread writes a unique filename so no lock is needed for the I/O.
    checkpoint_bytes = if checkpoint_dir !== nothing && !from_cache
        ShroffCelegansModels.write_dataset_checkpoint(checkpoint_dir, dataset.path, positions)
    else
        0
    end

    if from_cache && !count_done
        @debug "step6 dataset cache hit (already primed)" path=dataset.path group_idx=log_idx
    else
        n_done = Threads.atomic_add!(done_counter, 1) + 1
        @info("step6 dataset done",
            dataset = string(n_done, "/", grand_total),
            group_idx = log_idx,
            group_total = log_total,
            path = dataset.path,
            thread = Threads.threadid(),
            n_annotations = length(annotation_dict),
            n_timepoints = length(normalized_timepoints),
            elapsed_s = round(time() - ds_t0; digits=2),
            from_cache = from_cache,
            checkpoint_bytes = checkpoint_bytes,
        )
    end
    return positions
end

# Pin BLAS to 1 thread for a Threads.@threads region and restore afterward. The
# loop already parallelizes across all Julia threads, and each tps_solve! calls BLAS
# (bunchkaufman); if BLAS also multithreads we get nthreads × blas_threads contending
# for the same cores (on a 16-core pod: 16 Julia threads × a BLAS default of 50 ≈ 800
# threads on a 16-core quota). With nthreads()==1 there is no Julia parallelism, so
# leave BLAS as-is — a serial run still benefits from BLAS threads.
function _scope_blas_for_threaded_region()
    prev = LinearAlgebra.BLAS.get_num_threads()
    region = Threads.nthreads() > 1 ? 1 : prev
    LinearAlgebra.BLAS.set_num_threads(region)
    return prev, region
end

# Flat, all-datasets-at-once priming pass. Today the averaging loop parallelizes only
# over the ~3 datasets within one group (groups run sequentially), leaving a 16-core
# pod ~half idle. The per-dataset warp is independent across ALL datasets and ALL
# groups, and its output is the path-keyed `cache`, which the group pass already
# short-circuits on. So flatten every group's datasets into one Threads.@threads loop
# (16-wide) and fill the cache; the subsequent per-group averaging then runs entirely
# on cache hits. No change to numerical output — the cache value IS the computed value.
function prime_dataset_positions!(
    datasets::Dict{String, Vector{ShroffCelegansModels.Datasets.NormalizedDataset}},
    cache::Dict{String, Vector{Vector{Point3{Float64}}}},
    normalized_timepoints::AbstractVector{Float64} = LinRange(0,1,201);
    avg_models::Vector{<: CelegansModel} = avg_models,
    checkpoint_dir::Union{Nothing, AbstractString} = nothing,
)::Int
    @assert length(avg_models) == length(normalized_timepoints)

    # Flatten + dedupe by dataset.path (paths are group-independent; a duplicate
    # would be identical, and the cache is keyed by path anyway).
    seen = Set{String}()
    all_datasets = ShroffCelegansModels.Datasets.NormalizedDataset[]
    for k in keys(datasets), ds in datasets[k]
        if !(ds.path in seen)
            push!(seen, ds.path)
            push!(all_datasets, ds)
        end
    end

    # Model construction stays sequential (not run concurrently); only the warp is
    # threaded. Time it so the serial head is visible in the logs.
    @info "step6 prime: building datasets_info" n_datasets=length(all_datasets)
    t_info = time()
    datasets_info = get_datasets_info(all_datasets)
    @info "step6 prime: datasets_info done" elapsed_s=round(time() - t_info; digits=2)

    grand_total = length(datasets_info)
    done_counter = Threads.Atomic{Int}(0)
    cache_lock = ReentrantLock()
    prog = ProgressMeter.Progress(grand_total; desc="Prime annotations / dataset...")

    prev_blas, region_blas = _scope_blas_for_threaded_region()
    @info("step6 prime: BLAS threads scoped for @threads region",
        julia_nthreads = Threads.nthreads(),
        blas_in_region = region_blas,
        blas_prev = prev_blas,
    )
    try
        Threads.@threads for ds_idx in eachindex(datasets_info)
            compute_dataset_positions!(
                datasets_info[ds_idx], cache, cache_lock, normalized_timepoints, avg_models, checkpoint_dir;
                done_counter, grand_total,
                log_idx = ds_idx, log_total = grand_total,
                count_done = true,    # the prime pass is the authoritative counter
            )
            ProgressMeter.next!(prog)
        end
        ProgressMeter.finish!(prog)
    finally
        LinearAlgebra.BLAS.set_num_threads(prev_blas)
    end
    return grand_total
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

    # Global progress numbering: when the caller (the Dict-keyed average_annotations)
    # threads a shared counter + grand total through, the cache has already been
    # primed and the caller owns the count — so stay silent (@debug) on cache hits
    # here to avoid double-counting (count_done=false). Standalone calls (no counter)
    # fall back to this group's own dataset count and log every dataset.
    grand_total = progress_total === nothing ? length(datasets_info) : progress_total
    done_counter = progress_counter === nothing ? Threads.Atomic{Int}(0) : progress_counter
    count_done = progress_counter === nothing

    out = Vector{Vector{Dict{String, Point3{Float64}}}}(undef, length(datasets_info))
    cache_lock = ReentrantLock()
    prog = ProgressMeter.Progress(length(datasets_info); desc="Avg annotations / dataset...")
    prev_blas, region_blas = _scope_blas_for_threaded_region()
    @info("step6 BLAS threads scoped for @threads region",
        julia_nthreads = Threads.nthreads(),
        blas_in_region = region_blas,
        blas_prev = prev_blas,
    )
    try
    Threads.@threads for ds_idx in eachindex(datasets_info)
        dataset_info = datasets_info[ds_idx]
        positions = compute_dataset_positions!(
            dataset_info, cache, cache_lock, normalized_timepoints, avg_models, checkpoint_dir;
            done_counter, grand_total,
            log_idx = ds_idx, log_total = length(datasets_info),
            count_done = count_done,
        )
        annotation_dict = dataset_info.annotation_dict
        out[ds_idx] = try
            map(positions) do points
                Dict{String, Point3d}(keys(annotation_dict) .=> points)
            end
        catch err
            @error "Error creating dict for dataset $(dataset_info.dataset.path)" exception = (err, Base.catch_backtrace())
            map(positions) do _points
                Dict{String, Point3d}(keys(annotation_dict) .=> fill(Point3(NaN), length(keys(annotation_dict))))
            end
        end
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)
    finally
        LinearAlgebra.BLAS.set_num_threads(prev_blas)
    end
    return out
end
