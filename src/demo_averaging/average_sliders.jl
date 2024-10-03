function average_sliders(models::Vector{<:AbstractCelegansModel})
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    sliders = SliderGrid(f[2,1],
        (label="Weight", range=0:0.01:1),
        (label="Upsampling", range=0:4)
    )

    model = ShroffCelegansModels.average(models, AnalyticWeights([1.0, 0.0]), n_upsample=0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))

    on(throttle(0.1, sliders.sliders[1].value)) do value
        model = ShroffCelegansModels.average(models, AnalyticWeights([1-value, value]); n_upsample=sliders.sliders[2].value[])
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = ShroffCelegansModels.average(models, AnalyticWeights([1-value, value]); n_upsample=upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end

    display(f)
end

function average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    r = range(smts.modelTimeSeries.dataset.cell_key)
    sliders = SliderGrid(f[2,1],
        (label="Time (frames)", range=1.0:0.1:length(r)),
        (label="Upsampling", range=0:4)
    )

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

function record_average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries, movieFileName)
    f = average_sliders(smts)
    sg = f.content[2]
    set_close_to!(sg.sliders[2], 2)
    record(f, movieFileName, sg.sliders[1].range[]) do i
        set_close_to!(sg.sliders[1], i)
    end
end