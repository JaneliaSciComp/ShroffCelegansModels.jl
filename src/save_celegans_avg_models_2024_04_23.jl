h5open("celegans_avg_models_2024_04_23.h5", "w") do h5f
    for (i, m) in pairs(avg_models)
        save_celegans_model(h5f, @sprintf("avg_model_%03d", i), m)
    end
end

h5open("celegans_avg_models_2024_04_23.h5", "r+") do h5f
    attrs(h5f)["voxel_pitch_micrometers"] = 0.165
    h5f["measurements/volume"] = ShroffCelegansModels.volume_by_cross_section.(avg_models) .* 0.165^3
    get_length(model) = ShroffCelegansModels.central_spline(model)(1.0)[3]
    h5f["measurements/length"] = get_length.(avg_models) .* 0.165
    attrs(h5f["measurements/volume"])["units"] = "micrometers^3"
    attrs(h5f["measurements/length"])["units"] = "micrometers"
end