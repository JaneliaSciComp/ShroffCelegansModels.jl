function debug_annotation_ap_axis(
    avg_models,
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    second(x) = x[2]

    f = Figure(size = (1920, 1080))
    
    title_twisted = Observable("Twisted")
    lscene = LScene(f[1:3,1])
    ax_twisted = lscene.scene

    selected_annotation_name = Observable("")
    ax_distance = Axis(f[1,2],
        xlabel = "Arc distance along central spline",
        ylabel = "Distance to annotation",
        title = selected_annotation_name,
        limits = ((0,200), nothing)
    )
    ax_ratio = Axis(f[2,2],
        xlabel = "Arc distance along central spline",
        ylabel = "Ratio of distance to contour distance",
        title = selected_annotation_name,
        limits = ((0,200), nothing)
    )
    ax_z = Axis(f[3,2],
        xlabel = "Time (Normalized from Twitch to Hatch)",
        ylabel = "Z",
        title = "Z position of annotation",
        limits = ((0.0,1.0), nothing)
    )

    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[4, 1:2],
        (label="Time (Normalized)", range=r),
        (label="Exp. Factor", range=1.0:0.01:4),
    )
    slider_range = sliders.sliders[1].range[]

    annotation_text_toggle = Toggle(f)
    central_spline_lines_toggle = Toggle(f)
    contour_mesh_toggle = Toggle(f)
    f[5, 1:2] = grid!(
        [
            Label(f, "Contour mesh");;
            contour_mesh_toggle;;
            Label(f, "Annnotation text");;
            annotation_text_toggle;;
            Label(f, "Central spline lines");;
            central_spline_lines_toggle
        ]
    )

    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
        end
    end

    mts = smts.modelTimeSeries
    mts_nt = let _range = range(dataset.cell_key)
        x -> begin
            nt = x * (length(_range) - 1) + 1.0
            nt = round(Int, nt)
            title_twisted[] = "Twisted; idx = $nt, tp = $(_range[nt])"
            mts(nt)
        end
    end


    annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))
    twisted_central_spline = Observable(swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel)))))


    #=
    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = avg_models[idx]
        @info "annotation positions" nt
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end
    =#

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            if ismissing(ann)
                Point3f(NaN, NaN, NaN)
            else
                ann(nt)
            end
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    #=
    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end
    =#

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    #_annotation_cells = Observable(_annotation_positions_over_time[1])
    straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    twisted_positions = twisted_annotation_positions(0.0)
    twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    # @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    expansion_factor = sliders.sliders[2].value

    twisted_central_pts = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
    twisted_central_pts = Observable(twisted_central_pts)
    # @info twisted_central_pts[]

    twisted_central_line_match = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
    twisted_central_line_match = Observable(twisted_central_line_match)

    function get_annotation_text(dict)
        String.(get.((dataset.cell_key.mapping,), Symbol.(keys(dict)), String.(keys(dict))))
    end

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    #_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    #straight_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    #twisted_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions))))
    _annotation_text = get_annotation_text(annotation_dict)
    straight_annotation_text = get_annotation_text(annotation_dict)
    twisted_annotation_text = get_annotation_text(twisted_positions)
    twisted_annotation_text = Observable(twisted_annotation_text)

    annotation_menu = Menu(f[6, 1:2], options = twisted_annotation_text)

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    straight_color = Observable(repeat(color, length(smodel)))
    twisted_color = Observable(repeat(color, length(tmodel)))

    twisted_seam_cell_text = Observable(String.([tmodel.names[2:2:end]; tmodel.names[1:2:end]]))

    twisted_mesh_plot = mesh!(ax_twisted, twisted_mesh; colorrange, color = twisted_color, transparency = true, alpha = 0.5)
    connect!(twisted_mesh_plot.visible, contour_mesh_toggle.active)
    lines!(ax_twisted, twisted_central_spline)
    scatter!(ax_twisted, twisted_central_pts)
    cs_lines = lines!(ax_twisted, twisted_central_line_match)
    connect!(cs_lines.visible, central_spline_lines_toggle.active)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    ms_annotation_cells = meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    distances = Observable(Float64[])
    selected_distance = Observable(0.0)
    distance_central_pts = Observable(Point3f[])
    central_spline_arc_lengths = Observable(Float64[])

    selected_annotation_idx = Observable(1)
    max_r = ShroffCelegansModels.max_radius_function(tmodel)
    Npts = length(tmodel)
    z = LinRange(0, 1, Npts)
    # max_distance = Observable(max_r.(z) .* voxel_size .* first(sliders.sliders[2].range[]))
    max_distance = Observable(Float64[])
    lines!(ax_distance, central_spline_arc_lengths, distances)
    hlines!(ax_distance, selected_distance)
    lines!(ax_distance, central_spline_arc_lengths, max_distance)
    scatter!(ax_distance, central_spline_arc_lengths, distances, color = distances, colormap = Reverse(:viridis))
    scatter!(ax_twisted, distance_central_pts, color = distances, colormap = Reverse(:viridis))
    selected_twisted_annotation_cell = Observable([twisted_annotation_cells[][selected_annotation_idx[]]])
    meshscatter!(ax_twisted, selected_twisted_annotation_cell, color = :red, markersize=1.1)

    ratio = @lift try
        $distances ./ $max_distance
    catch err
        ones(size($max_distance))
    end
    lines!(ax_ratio, central_spline_arc_lengths, ratio)
    hlines!(ax_ratio, 1.0, linestyle = :dash)

    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        # title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        #_annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        # title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]

        tmodel = mts_nt(value)
        if !ismissing(tmodel)
            twisted_mesh[] = ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale)
            twisted_central_spline[] = swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel))))
            twisted_color[] = repeat(color, length(tmodel))
            twisted_seam_cells[] = swapyz_scale.(seam_cell_pts(tmodel, 0))
            twisted_seam_cell_labels[] = twisted_seam_cells[] .- Ref(Point3f(2,0,0))
            twisted_seam_cell_text[] = [String.(tmodel.names[1:2:end]); String.(tmodel.names[2:2:end])]
            annotation_positions = twisted_annotation_positions(value)
            if !ismissing(annotation_positions)
                twisted_annotation_text[] = get_annotation_text(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
            end
            #max_r = ShroffCelegansModels.max_radius_function(tmodel)
            #Npts = length(tmodel)
            #z = LinRange(0, 1, Npts)
            #max_distance[] = max_r.(z) .* voxel_size .* sliders.sliders[2].value[]
            #plot_distance(1)
            notify(annotation_menu.selection)
        end
    end

    on(throttle(0.1, sliders.sliders[2].value)) do expansion_factor_value
        value = sliders.sliders[1].value[]
        tmodel = mts_nt(value)
        annotation_positions = twisted_annotation_positions(value)
        if !ismissing(tmodel)
            if !ismissing(annotation_positions)
                twisted_annotation_text[] = get_annotation_text(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor_value)))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
                #max_r = ShroffCelegansModels.max_radius_function(tmodel)
                #Npts = length(tmodel)
                #z = LinRange(0, 1, Npts)
                #max_distance[] = max_r.(z) .* voxel_size .* expansion_factor_value
                plot_distance(selected_annotation_idx[])
            end
        end
    end



    function plot_distance(idx)
        selected_annotation_idx[] = idx
        #println(idx)
        selected_twisted_annotation_cell[] = [twisted_annotation_cells[][selected_annotation_idx[]]]
        cs = ShroffCelegansModels.central_spline(tmodel)
        Npts = length(tmodel)
        z = LinRange(0, 1, Npts)
        central_pts = swapyz_scale.(cs.(z))
        central_spline_arc_lengths.val = [0; cumsum(norm.(diff(central_pts)))]
        distance_central_pts[] = central_pts

        max_r = ShroffCelegansModels.max_radius_function(tmodel)
        expansion_factor_value = sliders.sliders[2].value[]
        max_distance[] = max_r.(z) .* voxel_size .* expansion_factor_value

        pt = twisted_annotation_cells[][idx]
        distances[] = norm.(central_pts .- pt)
        #autolimits!(ax_distance)
        #ylims!(ax_distance, nothing)
        limits!(ax_distance, (0, 200), (0, maximum(distances[])))
        selected_distance[] = norm(twisted_central_pts[][idx] - pt)

        selected_annotation_name[] = twisted_annotation_text[][idx]
    end
    plot_distance(::Nothing) = plot_distance(1)

    plot_distance(1)

    on(events(f).mousebutton, priority=2) do event
        # print("Mouse clicked")
        if event.button == Mouse.left && event.action == Mouse.press
            p, idx = pick(f)
            if p == ms_annotation_cells
                #println(twisted_annotation_text[][idx])
                plot_distance(idx)
            end
        end
    end

    z_positions = let idx=1
        Observable((x->x[idx][2]).(straight_annotation_positions_over_time))
    end
    lines!(ax_z, r, z_positions)
    DataInspector(ax_z)

    #common_annotations_text = collect(keys(annotation_dict))
    common_annotations_text = _annotation_text
    @info common_annotations_text

    on(annotation_menu.selection) do selected
        idx = findfirst(==(selected), twisted_annotation_text[])
        #if selected_annotation_idx[] != idx
            plot_distance(idx)
        #end
        if !isnothing(selected)
            idx_common = findfirst(==(selected), common_annotations_text)
            z_positions[] = (x->x[idx_common][2]).(straight_annotation_positions_over_time)
            autolimits!(ax_z)
        end
    end

    # display(f)
    print(typeof(straight_annotation_positions_over_time))
    f
end