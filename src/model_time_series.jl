struct ModelTimeSeries{Model <: AbstractCelegansModel, Cache <: LRU{<:AbstractString, Union{Model,Missing}}}
    dataset::NormalizedDataset
    cache::Cache
end
ModelTimeSeries(dataset::NormalizedDataset; max_cache_size = 2) =
    ModelTimeSeries(dataset, LRU{String, Union{CelegansModel,Missing}}(; maxsize = max_cache_size))
function (mts::ModelTimeSeries)(time_offset::Integer)
    try
        timepoint = range(mts.dataset.cell_key)[time_offset]
        lattice_filepath = get_lattice_filepath(mts.dataset, timepoint)
        ds = mts.dataset
        get!(mts.cache, lattice_filepath) do
            if timepoint ∈ ds.cell_key.outliers
                return missing
            elseif !isfile(lattice_filepath)
                @warn("$filepath is not a file on disk and is not marked as an outlier.")
                return missing
            else
                return build_celegans_model(lattice_filepath)
            end
        end
    catch err
        if err isa BoundsError
            return nothing
        else
            rethrow(err)
        end
    end
end

struct StraightenedModelTimeSeries{Model <: AbstractCelegansModel, Cache <: LRU{<: AbstractString, Union{Model,Missing}}, ModelTS <: ModelTimeSeries}
    modelTimeSeries::ModelTS
    cache::Cache
end
function StraightenedModelTimeSeries(dataset::NormalizedDataset, model_cache::LRU, straightened_model_cache::LRU)
    mts = ModelTimeSeries(dataset, model_cache)
    return StraightenedModelTimeSeries(mts, straightened_model_cache)
end
function StraightenedModelTimeSeries(dataset::NormalizedDataset; max_cache_size = 2)
    mts = ModelTimeSeries(dataset; max_cache_size)
    StraightenedModelTimeSeries(mts, LRU{String, Union{StraightenedCelegansModel,Missing}}(; maxsize = max_cache_size))
end
function (smts::StraightenedModelTimeSeries)(time_offset::Integer)
    mts = smts.modelTimeSeries
    try
        timepoint = range(mts.dataset.cell_key)[time_offset]
        lattice_filepath = get_lattice_filepath(mts.dataset, timepoint)
        get!(smts.cache, lattice_filepath) do
            model = mts(time_offset)
            if ismissing(model)
                missing
            else
                StraightenedCelegansModel(model)
            end
        end
    catch err
        if err isa BoundsError
            return nothing
        else
            rethrow(err)
        end
    end
end
function (smts::StraightenedModelTimeSeries)(time_offset; n_upsample = 0)
    next = nextModelIndex(smts, time_offset)
    prev = prevModelIndex(smts, time_offset)
    if isnothing(next) && isnothing(prev)
        return nothing
    end
    if isnothing(next)
        return average([smts(prev)]; n_upsample)
    end
    if isnothing(prev)
        return average([smts(next)]; n_upsample)
    end
    if next == prev
        return average([smts(next)]; n_upsample)
    end
    return average([smts(prev), smts(next)], AnalyticWeights([next - time_offset, time_offset - prev]); n_upsample)
end

function get_lattice(ds::NormalizedDataset, time_offset=1)::Union{Missing, String}
    timepoint = range(ds.cell_key)[time_offset]
    if timepoint ∈ ds.cell_key.outliers
        return missing
    else
        filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
        if isfile(filepath)
            return filepath
        else
            @warn("$filepath is not a file on disk and is not marked as an outlier.")
            return missing
        end
    end
end

function get_lattice_filepath(ds::NormalizedDataset, timepoint::Integer)::String
    filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
    return filepath
end

function nextModelIndex(smts::StraightenedModelTimeSeries, index)
    index = ceil(Int, index)
    while ismissing(smts(index))
        index += 1
    end
    return index
end

function prevModelIndex(smts::StraightenedModelTimeSeries, index)
    index = floor(Int, index)
    while ismissing(smts(index))
        index -= 1
    end
    return index
end