#const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_full.json"
if gethostname() == "KITTISOPIKULM-2"
    const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_2026_03_19_v2.json"
else
    const config_path = joinpath(@__DIR__, "..", "..", "config", "linux", "config_2026_03_19_v2.json")
end
const voxel_size = 0.1625 # um

using LinearAlgebra
using ShroffCelegansModels.JSON3
LinearAlgebra.BLAS.set_num_threads(12)

@info "Reading Config JSON"
config_json, cell_keys, datasets = read_config_json()

flattened_datasets = collect(Iterators.flatten(values(datasets)))

min_length = map(flattened_datasets) do ds
    length(range(ds.cell_key))
end |> minimum

# models, smodels = build_models_over_time(datasets["RW10598"][2]);

normal_flattened_datasets = filter(x->length(range(x.cell_key)) != 259, flattened_datasets)
normal_flattened_datasets = filter(x->x.cell_key.name != "Vab-1_Pos0",normal_flattened_datasets)
smts_datasets = ShroffCelegansModels.StraightenedModelTimeSeries.(normal_flattened_datasets)
lengths = normal_flattened_datasets .|> x->length(range(x.cell_key))
smts_datasets_nt = map(zip(smts_datasets, lengths)) do (ds, _length)
    x -> begin
        nt = x * (_length - 1) + 1.0
        ds(nt)
    end
end
models_at_nt(nt) = map(smts_datasets_nt) do ds
    ds(nt)
end
include("../save_celegans_avg_models.jl")

recalculate_avg_models = false

if recalculate_avg_models
    @info "Calculating average models"
    avg_models = get_avg_models()
else
    @info "Loading average models"
    avg_models = load_avg_models("celegans_avg_models_2024_07_26.h5")
end

using ShroffCelegansModels: annotations_cache, annotation_position_cache, my_annotation_position_cache

    int_ds = filter(flattened_datasets) do ds
        "int1dr" in values(ds.cell_key.mapping)
    end
    smts2 = ShroffCelegansModels.StraightenedModelTimeSeries(int_ds[2])
    _length2 = length(range(int_ds[2].cell_key))
    smts_nt2 = x -> begin
        nt = x * (_length2 - 1) + 1.0
        smts2(nt, 2)
    end
