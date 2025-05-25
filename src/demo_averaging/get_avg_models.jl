function get_avg_models(n=201)
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