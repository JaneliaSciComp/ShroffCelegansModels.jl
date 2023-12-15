using ShroffCelegansModels
lattice_final_path = joinpath(pkgdir(ShroffCelegansModels), "test", "artifacts", "lattice.csv")
        model = ShroffCelegansModels.build_celegans_model(lattice_final_path)
        smodel = ShroffCelegansModels.StraightenedCelegansModel(model)
        get_model_contour_mesh(model)
        get_model_contour_mesh(smodel)