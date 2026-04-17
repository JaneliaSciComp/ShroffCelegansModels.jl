using ShroffCelegansModels
using ProgressMeter

function get_avg_models(n=201)
    config_json, cell_keys, datasets = read_config_json()

    flattened_datasets = collect(Iterators.flatten(values(datasets)))

    min_length = map(flattened_datasets) do ds
        length(range(ds.cell_key))
    end |> minimum


    normal_flattened_datasets = filter(x->length(range(x.cell_key)) != 259, flattened_datasets)
    normal_flattened_datasets = filter(x->x.cell_key.name != "Vab-1_Pos0",normal_flattened_datasets)

    smts_datasets = ShroffCelegansModels.StraightenedModelTimeSeries.(normal_flattened_datasets)
    lengths = normal_flattened_datasets .|> x->length(range(x.cell_key))
    smts_datasets_nt = map(zip(smts_datasets, lengths)) do (ds, _length)
        x -> begin
            nt = x * (_length - 1) + 1.0
            # @info "Normalized time" nt
            ds(nt)
        end
    end
    models_at_nt(nt) = map(smts_datasets_nt) do ds
        ds(nt)
    end


    r = LinRange(0.0, 1.0, n)
    first_avg_model = let models = models_at_nt(r[1])
        models = filter(!isnothing, models)
        models = identity.(models)
        ShroffCelegansModels.average(models; n_upsample = 2)
    end
    avg_models = Vector{typeof(first_avg_model)}(undef, length(r))
    avg_models[1] = first_avg_model
    @showprogress desc="Averaging models..." Threads.@threads for i in eachindex(r)[2:end]
        nt = r[i]
        models = models_at_nt(nt) 
        models = filter(!isnothing, models)
        models = identity.(models)
        avg_models[i] = ShroffCelegansModels.average(models; n_upsample = 2)
    end
    return avg_models

    #=
    avg_models = map(r) do nt
        @info nt
        models = models_at_nt(nt)
        models = filter(!isnothing, models)
        models = identity.(models)
        ShroffCelegansModels.average(models; n_upsample = 2)
    end
    =#
end