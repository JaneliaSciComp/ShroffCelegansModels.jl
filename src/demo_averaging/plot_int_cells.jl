const int_cells = ["Mu_int_L", "int1dl", "int1dr", "int1vl", "int1vr", "int2d", "int2v", "int3d", "int3v", "int4d", "int4v", "int5l", "int5r", "int6l", "int6r", "int7l", "int7r", "int8l", "int8r", "int9l", "int9r"]

function plot_int_cells(avg_models; flattened_datasets = flattened_datasets)
    int_ds = filter(flattened_datasets) do ds
        "int1dr" in values(ds.cell_key.mapping)
    end

    # models, smodels = build_models_over_time(int_ds[1]);
    
    ds_idx = 1

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(int_ds[ds_idx])
    _length = length(range(int_ds[ds_idx].cell_key))
    smts_nt = x -> begin
        nt = x * (_length - 1) + 1.0
        # @info "Normalized time" nt
        smts(nt, 2)
    end

    smodel = smts_nt(1.0)

    timepoint = 87
    int_annotations_1 = load_straightened_annotations_over_time(int_ds[ds_idx])

    warp_from = ShroffCelegansModels.lattice(smodel) |> vec
    warp_to = ShroffCelegansModels.lattice(avg_models[end]) |> vec
    @info "types" typeof(warp_from) typeof(warp_to)
    tps_solved = tps_solve(warp_from, warp_to, 1)

    pts = collect(values(int_annotations_1[timepoint]))
    @info "pts" typeof(pts)
    transformed_pts = ThinPlateSplines.tps_deform(pts, tps_solved)

    pts = swapyz_scale.(pts)
    transformed_pts = swapyz_scale.(transformed_pts)

    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)

    mesh!(ax, smodel; transform_points = swapyz_scale, alpha = 0.1, transparency = true, color = :green)
    mesh!(ShroffCelegansModels.get_model_contour_mesh(avg_models[timepoint], transform_points = swapyz_scale), alpha = 0.1, transparency = true, color = :magenta)

    scatter!(ax, pts, color = :green)
    scatter!(ax, transformed_pts, color = :magenta)
    
    text!(ax, pts; text = collect(keys(int_annotations_1[timepoint])))
    text!(ax, transformed_pts; text = collect(keys(int_annotations_1[timepoint])))

    arrows!(ax, pts, transformed_pts .- pts)

    ylims!(ax, (0, 1000*voxel_size))
    return f
end

