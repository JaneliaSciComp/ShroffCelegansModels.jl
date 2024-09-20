function debug_average_models_with_annotations(
    avg_models,
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    second(x) = x[2]

    f = Figure(size = (1920, 1080))
    
    title_prewarp = Observable("Prewarp")
    ax_prewarp = Axis3(f[1,2], aspect = (1, 10, 1); title = title_prewarp)
    
    title = Observable("Title")
    ax = Axis3(f[2,2], aspect = (1, 10, 1); title)

    title_twisted = Observable("Twisted")
    ax_twisted = Axis3(f[1:2,1]; title = title_twisted)

    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[3,1:2],
        (label="Time (Normalized)", range=r),
        (label="Exp. Factor", range=1.0:0.01:4),
    )
    slider_range = sliders.sliders[1].range[]

    annotation_text_toggle = Toggle(f)
    central_spline_lines_toggle = Toggle(f)
    f[4, 1:2] = grid!(
        [Label(f, "Annnotation text") annotation_text_toggle Label(f, "Central spline lines") central_spline_lines_toggle]
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

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#

    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))
    twisted_central_spline = Observable(swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel)))))


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

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            ann(nt)
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    _annotation_cells = Observable(_annotation_positions_over_time[1])
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


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    straight_color = Observable(repeat(color, length(smodel)))
    twisted_color = Observable(repeat(color, length(tmodel)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    ann_txt = text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))

    mesh!(ax_prewarp, straight_mesh; colorrange, color = straight_color, transparency = true, alpha = 0.1)
    meshscatter!(ax_prewarp, straight_seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax_prewarp, straight_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax_prewarp, straight_seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    ann_txt = text!(ax_prewarp, straight_annotation_cells; text = straight_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    #lines!(ax_prewarp, _lines, color = :black)
    ylims!(ax_prewarp, (0, 200))

    twisted_seam_cell_text = Observable(String.([tmodel.names[2:2:end]; tmodel.names[1:2:end]]))

    mesh!(ax_twisted, twisted_mesh; colorrange, color = twisted_color, transparency = true, alpha = 0.5)
    lines!(ax_twisted, twisted_central_spline)
    scatter!(ax_twisted, twisted_central_pts)
    cs_lines = lines!(ax_twisted, twisted_central_line_match)
    connect!(cs_lines.visible, central_spline_lines_toggle.active)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        _annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
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
            end
        end
    end

    display(f)
    f
end

function debug_average_models_with_annotations(
    avg_models,
    datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset};
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    dataset = datasets[1]

    f = Figure(size = (1920, 1080))
    
    title_prewarp = Observable("Prewarp")
    ax_prewarp = Axis3(f[1,2], aspect = (1, 10, 1); title = title_prewarp)
    
    title = Observable("Title")
    ax = Axis3(f[2,2], aspect = (1, 10, 1); title)

    title_twisted = Observable("Twisted")
    ax_twisted = Axis3(f[1:2,1]; title = title_twisted)

    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[3,1:2],
        (label="Time (Normalized)", range=r),
        (label="Exp. Factor", range=1.0:0.01:1.5)
    )
    slider_range = sliders.sliders[1].range[]

    dataset_paths = map(x->x.path, datasets)
    dataset_menu = Menu(f[4,1:2], options = zip(dataset_paths, datasets), default = datasets[2].path)

    dataset = dataset_menu.selection[]

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

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#

    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))

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

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            ann(nt)
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    twisted_positions = twisted_annotation_positions(0.0)
    twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    _annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))
    straight_annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))
    twisted_annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions)))


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    straight_color = Observable(repeat(color, length(smodel)))
    twisted_color = Observable(repeat(color, length(tmodel)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))

    mesh!(ax_prewarp, straight_mesh; colorrange, color = straight_color, transparency = true, alpha = 0.1)
    meshscatter!(ax_prewarp, straight_seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax_prewarp, straight_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax_prewarp, straight_seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    text!(ax_prewarp, straight_annotation_cells; text = straight_annotation_text, align = (:right, :bottom))
    #lines!(ax_prewarp, _lines, color = :black)
    ylims!(ax_prewarp, (0, 200))

    twisted_seam_cell_text = Observable([tmodel.names[2:2:end]; tmodel.names[1:2:end]])

    mesh!(ax_twisted, twisted_mesh; colorrange, color = twisted_color, transparency = true, alpha = 0.5)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    on(dataset_menu.selection) do ds
        dataset = ds
        smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
        smts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                smts(nt, 2)
            end
        end

        mts = smts.modelTimeSeries
        mts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                nt = round(Int, nt)
                title_twisted[] = "Twisted; idx = $nt, tp = $(_range[nt])"
                mts(nt)
            end
        end


        annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

        twisted_annotation_positions = let _length = length(range(dataset.cell_key))
        function twisted_annotation_positions(nt)
            idx = round(Int, nt*(_length - 1) + 1)
            dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
            Dict(keys(dict) .=> swapyz_scale.(values(dict)))
        end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    end

    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        _annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]

        tmodel = mts_nt(value)
        if !ismissing(tmodel)
            twisted_mesh[] = ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale)
            twisted_color[] = repeat(color, length(tmodel))
            twisted_seam_cells[] = swapyz_scale.(seam_cell_pts(tmodel, 0))
            twisted_seam_cell_labels[] = twisted_seam_cells[] .- Ref(Point3f(2,0,0))
            twisted_seam_cell_text[] = [tmodel.names[1:2:end]; tmodel.names[2:2:end]]
            annotation_positions = twisted_annotation_positions(value)
            if !ismissing(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
            end
        end
    end

    display(f)
    f
end