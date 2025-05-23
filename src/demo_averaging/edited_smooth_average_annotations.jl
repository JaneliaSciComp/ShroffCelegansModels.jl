smoothed_average_annotations = smooth_average_annotations(
    average_annotations_dict;
    smooth_factor_r = 0.20,
    smooth_factor_θ = 0.20,
    smooth_factor_z = 0.30
)
f = isnothing.(match.(r"P\d/\d*[LR]", smoothed_average_annotations["RW10742"].annotations))
edited_smoothed_average_annotations = copy(smoothed_average_annotations)
edited_smoothed_average_annotations["RW10742"] = (;
    annotations=edited_smoothed_average_annotations["RW10742"].annotations[f],
    positions=map(smoothed_average_annotations["RW10742"].positions) do p
        p[f]
    end
)
delete!(edited_smoothed_average_annotations, "RW10896")
save_average_annotations(edited_smoothed_average_annotations; filename = "edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells_2025_02_13.h5")

for (k, v) in edited_smoothed_average_annotations
    if "p9/10l" in lowercase.(v.annotations)
        println(k)
    end
end