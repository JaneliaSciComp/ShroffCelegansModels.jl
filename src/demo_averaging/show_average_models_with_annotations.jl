function show_average_models_with_annotations(
    models,
    dataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    #vol_ax = Axis(f[2,1]; xlabel = "Time (Normalized)", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    #length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[2,1],
        (label="Time (Normalized)", range=r),
    )
    slider_range = sliders.sliders[1].range[]

    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
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


    model = models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)

    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = models[idx]
        @info "annotation positions" nt
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end

    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    # _annotation_cells = Observable(annotation_positions(0.0))

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    _annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))


    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        #_lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        #title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        #current_volume[] = Point(value, volumes[idx])
        #current_length[] = Point(value, lengths[idx])
        #smodel = smts_nt(value)
        #=
        _annotation_cells[] = swapyz_scale.(transform_annotations(
            smodel, model, map(values(annotation_dict)) do ann
                ann(value)
            end
        ))
        =#
        _annotation_cells[] = _annotation_positions_over_time[idx]
    end

    display(f)
    f
end

function record_show_average_models_with_annotations(
    models,
    dataset,
    movieFileName;
    use_myuntwist = false
)
    f = show_average_models_with_annotations(models, dataset; use_myuntwist)
    sg = f.content[2]
    record(f, movieFileName, sg.sliders[1].range[]; framerate = 10) do i
        set_close_to!(sg.sliders[1], i)
    end
end
