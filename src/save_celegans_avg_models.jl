using Dates

avg_models_filename = "celegans_avg_models_" * Dates.format(now(), "yyyy_mm_dd") * ".h5"

function save_avg_models(avg_models_filename = avg_models_filename, avg_models = avg_models)
    h5open(avg_models_filename, "w") do h5f
        for (i, m) in pairs(avg_models)
            save_celegans_model(h5f, @sprintf("avg_model_%03d", i), m)
        end
    end
end

function load_avg_models(avg_models_filename= avg_models_filename)
    avg_models = []
    h5open(avg_models_filename, "r") do h5f
        for k in keys(h5f)
            group_attrs = attrs(h5f[k])
            if haskey(group_attrs, "julia_type") &&
               contains(group_attrs["julia_type"], "ShroffCelegansModels.Types.CelegansModel")
                push!(avg_models, load_celegans_model(h5f[k]))
            end
        end
    end
    identity.(avg_models)
end

function save_measurements(avg_models_filename = avg_models_filename)
    h5open(avg_models_filename) do h5f
        attrs(h5f)["voxel_pitch_micrometers"] = 0.165
        h5f["measurements/volume"] = ShroffCelegansModels.volume_by_cross_section.(avg_models) .* 0.165^3
        get_length(model) = ShroffCelegansModels.central_spline(model)(1.0)[3]
        h5f["measurements/length"] = get_length.(avg_models) .* 0.165
        attrs(h5f["measurements/volume"])["units"] = "micrometers^3"
        attrs(h5f["measurements/length"])["units"] = "micrometers"
    end
end