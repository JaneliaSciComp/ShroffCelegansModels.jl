function average_sliders(models::Vector{<:AbstractCelegansModel})
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    sliders = SliderGrid(f[2,1],
        (label="Weight", range=0:0.01:1, visible=false),
        (label="Upsampling", range=0:4, visible=false)
    )

    model = ShroffCelegansModels.average(models, AnalyticWeights([1.0, 0.0]), n_upsample=0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    _lines = Observable(swapyz_scale.(cross_sections_at_knots(model)))

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
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _lines[] = swapyz_scale.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = ShroffCelegansModels.average(models, AnalyticWeights([1-value, value]); n_upsample=upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _lines[] = swapyz_scale.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end

    display(f)
end

function average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    #ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    lscene = LScene(f[1,1], show_axis=false) 
    ax = lscene.scene
    r = range(smts.modelTimeSeries.dataset.cell_key)
    sliders = SliderGrid(f[2,1],
        (label="Time (frames)", range=1.0:0.1:length(r)),
        (label="Upsampling", range=0:4)
    )

    # Time Label
    label_offset = 12 
    time_text = Observable("t = 14:00")
    text!(-label_offset, 0, label_offset; text = time_text)

    # Scalebar
    scalebar_size_um = 10
    scalebar_y_offset = 160
    scalebar_y_offset_obs = Observable(180)
    scalebar_text = Observable("$scalebar_size_um μm")
    scalebar = Observable(Point3f[
        [label_offset+1, scalebar_y_offset,                    -label_offset-1],
        [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1]
    ])
    scalebar_label_position = lift(scalebar) do pos
        first(pos) + Point3f(0,0,1)
    end
    text!( scalebar_label_position; text = scalebar_text)
    lines!(scalebar, color = :white, linewidth = 5)

    model = smts(length(r), 0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    _lines = Observable(swapyz_scale.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, 0)))
    # _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))
    # _seam_cell_labels = Observable(_seam_cells[] .+ Ref(Point3f(10,0,0)))
    _seam_cell_labels = _seam_cells
    right_seam_cell_labels = lift(_seam_cells) do sc
        sc[1:end÷2] .+ Point3f(10,1,2)
    end
    left_seam_cell_labels = lift(_seam_cells) do sc
        sc[end÷2+1:end] .+ Point3f(-10,1,2)
    end

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading, visible=true)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 0.5)
    #text!(ax, _seam_cell_labels; text = [replace.(model.names[1:2:end], 'L' => 'R'); model.names[1:2:end]], align = (:right, :bottom))
    text!(ax, right_seam_cell_labels;
        text = replace.(model.names[1:2:end], 'L' => 'R'),
        align = (:left, :bottom),
    )
    text!(ax, left_seam_cell_labels;
        text = model.names[1:2:end],
        align = (:left, :bottom),
        visible=false
    )
    lines!(ax, _lines, color = :black, visible=false)
    #ylims!(ax, (0, 1000))

    #on(throttle(0.1, sliders.sliders[1].value)) do value
    on(sliders.sliders[1].value) do value
        total_minutes = (value-1)*5
        println(total_minutes)
        hours = round(Int, total_minutes/60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        time_text[] = "t = $hours:$(@sprintf("%02d", minutes))"
        n_upsample = sliders.sliders[2].value[]
        model = smts(value, n_upsample)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _lines[] = swapyz_scale.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        #_seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = smts(value, upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _lines[] = swapyz_scale.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, upsampling))
        # _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
        # _seam_cell_labels[] = _seam_cells[] .+ Ref(Point3f(10,0,0))
    end

    display(f)

    cc = Camera3D(ax;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(60, 90, 0)
    )

    #zoom!(ax, cc, 0.5)

    f
end

function record_average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries, movieFileName; interpolate = false, kwargs...)
    f = average_sliders(smts)
    sg = f.content[2]
    # hide slider grid
    for child in sg.blockscene.children
           child.visible = false
    end
    set_close_to!(sg.sliders[2], 2)
    _range = sg.sliders[1].range[]
    if !interpolate
        # Do not interpolate in this video
        _range = first(_range):last(_range)
    end
    record(f, movieFileName, _range; update=false, kwargs...) do i
        set_close_to!(sg.sliders[1], i)
        println(i)
    end
end