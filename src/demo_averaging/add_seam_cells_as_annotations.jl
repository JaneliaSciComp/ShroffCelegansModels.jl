function seam_cells_as_annotations(avg_models)
    seam_cell_names = [
        replace.(avg_models[1].names[1:2:end], "L" => "R");
        avg_models[1].names[1:2:end]
    ]
    return (;
        annotations = seam_cell_names,
        positions = map(avg_models) do avg_model
            swapyz_scale.(seam_cell_pts(avg_model, 2))
        end
    )
end
# average_annotations_dict["seam_cells"] = seam_cells_as_annotations(avg_models)