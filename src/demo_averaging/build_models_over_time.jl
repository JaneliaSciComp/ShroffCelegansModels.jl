function build_models_over_time(dataset::ShroffCelegansModels.Datasets.NormalizedDataset, offsets = 1:length(range(dataset.cell_key)))
    models = map(offsets) do time_offset
        lattice = get_lattice(dataset, time_offset)
        ShroffCelegansModels.build_celegans_model(lattice)
    end
    smodels = map(models) do model
        if ismissing(model)
            missing
        else
            ShroffCelegansModels.Types.StraightenedCelegansModel(model)
        end
    end
    return models, smodels
end