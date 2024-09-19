function read_config_json(config_path::AbstractString = config_path)
    config_json = JSON3.read(config_path)
    cell_keys = Dict{String, Vector{ShroffCelegansModels.CellKey}}()
    map(config_json.data.strains) do strain
        cell_keys[strain.name] = map(strain.folderpaths) do folder_path
            cell_key_path = joinpath(folder_path, "cell_key.json")
            try
                cell_key_json = JSON3.read(cell_key_path)
                ShroffCelegansModels.CellKey(cell_key_json)
            catch err
                if isfile(cell_key_path)
                    throw(ArgumentError("Problem parsing $cell_key_path"))
                else
                    throw(ErrorException("$cell_key_path does not exist"))
                end
            end
        end
    end

    datasets = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
    map(config_json.data.strains) do strain
        datasets[strain.name] = map(strain.folderpaths) do folder_path
            ShroffCelegansModels.NormalizedDataset(joinpath(folder_path, "RegB"))
        end
    end
    return config_json, cell_keys, datasets
end