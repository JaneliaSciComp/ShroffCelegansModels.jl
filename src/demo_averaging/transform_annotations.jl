function transform_annotations(from_model, to_model, annotations::AbstractVector{<: Point})
    warp_from = ShroffCelegansModels.lattice(from_model) |> vec
    warp_to = ShroffCelegansModels.lattice(to_model) |> vec
    tps_solved = tps_solve(warp_from, warp_to, 1)

    transformed_pts = ThinPlateSplines.tps_deform(annotations, tps_solved)
    return transformed_pts
end
function transform_annotations(from_model, to_model, annotations::Dict)
    return Dict(keys(annotations) .=> transform_annotations(from_model, to_model, collect(values(annotations))))
end