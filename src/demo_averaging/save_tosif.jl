# resave_for_tosif("edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5")
function resave_for_tosif(filename; target_filename = replace(filename, ".h5" => "_tosif.h5"))
    if isfile(target_filename)
        error("$target_filename exists")
    end
    tosif_data = h5open(filename) do h5f
        strain_matrices = Array{Float64,3}[]
        annotation_labels = String[]
        strain_labels = String[]
        for strain in keys(h5f)
            matrices = Matrix{Float64}[]
            for tp in 1:201
                tpk = @sprintf("timepoint_%03d", tp)
                ds = h5f[strain][tpk]
                push!(matrices, ds[])
            end
            labels = h5f[strain]["annotations"][]
            append!(annotation_labels, labels)
            append!(strain_labels, repeat([strain], length(labels)))
            push!(strain_matrices, stack(matrices, dims=2))
        end
        minutes = (0:0.005:1) * 420 .+ 420 |> collect
        trajectories = cat(strain_matrices..., dims=1)
        trajectories = trajectories[:,:,[1,3,2]]
        (; trajectories, annotation_labels, strain_labels, minutes)
    end

    h5open(target_filename, "w") do h5f
        h5f["annotations"] = tosif_data.annotation_labels
        h5f["strains"] = tosif_data.strain_labels
        h5f["time_minutes"] = tosif_data.minutes
        h5f["trajectories"] = tosif_data.trajectories
    end
    return tosif_data
end