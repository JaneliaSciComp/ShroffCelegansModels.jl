function average_sliders_with_volumes(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    vol_ax = Axis(f[2,1]; xlabel = "Frames", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    r = range(smts.modelTimeSeries.dataset.cell_key)
    sliders = SliderGrid(f[3,1],
        (label="Time (frames)", range=1.0:0.1:length(r)),
        (label="Upsampling", range=0:4)
    )
    slider_range = sliders.sliders[1].range[]

    # precalculate the volumes
    volumes_and_lengths = map(slider_range) do t
        model = smts(t)
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)

    # lengths
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)

    model = smts(1.0, 0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz.(seam_cell_pts(model, 0)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    meshscatter!(ax, _seam_cells; markersize = 10.0, color = :gray, alpha = 0.5)
    text!(ax, _seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))

    on(throttle(0.1, sliders.sliders[1].value)) do value
        n_upsample = sliders.sliders[2].value[]
        model = smts(value, n_upsample)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
        current_volume[] = Point(value, volumes[floor(Int, (value-1)/step(slider_range) + 1)])
        current_length[] = Point(value, lengths[floor(Int, (value-1)/step(slider_range) + 1)])
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = smts(value, upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, upsampling))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
    end

    display(f)
    f
end

function record_average_sliders_with_volumes(smts::ShroffCelegansModels.StraightenedModelTimeSeries, movieFileName)
    f = average_sliders_with_volumes(smts)
    sg = f.content[4]
    set_close_to!(sg.sliders[2], 2)
    record(f, movieFileName, sg.sliders[1].range[]) do i
        set_close_to!(sg.sliders[1], i)
    end
end