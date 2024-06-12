function get_lattice(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)::Union{Missing, String}
    timepoint = range(ds.cell_key)[time_offset]
    if timepoint ∈ ds.cell_key.outliers
        return missing
    else
        filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
        if isfile(filepath)
            return filepath
        else
            throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
        end
    end
end