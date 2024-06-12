#const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_full.json"
const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_2024_01_10_v2.json"
const voxel_size = 0.1625 # um

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(12)

config_json, cell_keys, datasets = read_config_json()

flattened_datasets = collect(Iterators.flatten(values(datasets)))

min_length = map(flattened_datasets) do ds
    length(range(ds.cell_key))
end |> minimum

models, smodels = build_models_over_time(datasets["RW10598"][2]);

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
# r = LinRange(0.0, 1.0, 201)
#=
avg_models = map(r) do nt
    @info nt
    models = models_at_nt(nt)
    models = filter(!isnothing, models)
    models = identity.(models)
    ShroffCelegansModels.average(models; n_upsample = 2)
end
=#

avg_models = get_avg_models()


#=
nothing_count = map(r) do nt
    @info nt
    models = models_at_nt(nt)
    return sum(isnothing.(models))
    models = filter(!isnothing, models)
    models = identity.(models)
    ShroffCelegansModels.average(models; n_upsample = 2)
end
=#


const annotation_position_cache = Dict{String, Any}()
const my_annotation_position_cache = Dict{String, Any}()

    int_ds = filter(flattened_datasets) do ds
        "int1dr" in values(ds.cell_key.mapping)
    end
    smts2 = ShroffCelegansModels.StraightenedModelTimeSeries(int_ds[2])
    _length2 = length(range(int_ds[2].cell_key))
    smts_nt2 = x -> begin
        nt = x * (_length2 - 1) + 1.0
        # @info "Normalized time" nt
        smts2(nt, 2)
    end

