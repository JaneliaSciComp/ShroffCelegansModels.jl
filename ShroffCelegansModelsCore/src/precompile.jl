@setup_workload begin
    lattice_final_path = joinpath(pkgdir(ShroffCelegansModelsCore), "test", "fixtures", "lattice.csv")
    @compile_workload begin
        try
            model = ShroffCelegansModelsCore.build_celegans_model(lattice_final_path)
            smodel = ShroffCelegansModelsCore.StraightenedCelegansModel(model)
            get_model_contour_mesh(model)
            get_model_contour_mesh(smodel)
        catch err
            Base.showerror(stdout, err, Base.catch_backtrace())
        end
    end
end