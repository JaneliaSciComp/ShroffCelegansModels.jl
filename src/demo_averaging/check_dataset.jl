function check_dataset(dataset; show_data_frames = false)
    r = range(dataset.cell_key)
    for t in r
        filepath = ShroffCelegansModels.get_lattice_filepath(dataset, t)
        if !isfile(filepath)
            if t ∉ dataset.cell_key.outliers
                @error "$filepath does not exist and is not an outlier"
                println()
            end
            continue
        end
        try
            df = CSV.read(filepath, DataFrame)
            _names = df[:,1]
            repeated_names = String[]
            if !allunique(df, 1)
                prev_name = ""
                for name in sort(_names)
                    if name == prev_name
                        push!(repeated_names, name)
                    end
                    prev_name = name
                end
            end
            missing_names = String[]
            for psc in ShroffCelegansModels.psuedo_seam_cells
                if psc ∉ _names
                    push!(missing_names, psc)
                end
                psc_R = replace(psc, 'L' => 'R')
                if psc_R ∉ _names
                    push!(missing_names, psc_R)
                end
            end
            if !isempty(repeated_names) || !isempty(missing_names)
                if !isempty(repeated_names) && !isempty(missing_names)
                    @error "Name error in $filepath" repeated_names missing_names
                elseif !isempty(repeated_names)
                    @error "Repeated names in $filepath" repeated_names
                elseif !isempty(missing_names)
                    @error "Missing names in $filepath" missing_names
                end
                show_data_frames && display(df)
                println()
            end
        catch err
            @error "Error parsing $filepath" err
        end
    end
end