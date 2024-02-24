using GLMakie

raw_colors = [
    (5, 141, 150)
    (0, 164, 80)
    (82, 180, 72)
    (138, 195, 65)
]
colors = map(raw_colors) do int_color
    RGB((int_color ./ 255)...)
end

fig_size = (1920, 1080)


f = Figure(
    size=fig_size,
    backgroundcolor = :lightgray
)

left_panel = f[1,1] = GridLayout()
center_panel = f[1,2] = GridLayout()
right_panel = f[1,3] = GridLayout()

Label(
    left_panel[1,1],
    text = "Abstract",
    width = 640,
    fontsize = 24
)

Label(
    left_panel[2,1],
    text =
"""
Caenorhabditis elegans, a roundworm, is a laboratory model organism that is amenable to the study of cellular and systems level physiology. By imaging the C. elegans embryo using dual view Selective Plane Illumination Microscopy (diSPIM) in three dimensions over time, we are developing a spatial temporal model of early nematode embryogenesis. A three dimensional lattice representation is generated using MIPAV (Medical Image Processing, Analysis and Visualization), a Java-based visualization program, by tracking the positions of seam cells from twitch to hatch in time-lapse volumetric images. We are adapting MIPAV to be able to refine the lattice by editing individual cross sections to better match the physical constraints of the worm. The result is a spatially continuous model interpolated using natural cubic B-splines and Fourier series. Furthermore, we are developing software in Julia to untwist the 3D lattice model of individual worms by transforming it into a cylindrical coordinate system relative to the organism's central anterior-posterior axis. These models are averaged to create a continuous meshwork over time through linear interpolation and registered together to create a common reference model. The relative positions of tracked, fluorescently labeled cells are then collected into a common 4D cellular atlas of C. elegans embryogenesis. Overall, we have developed the next generation of our C. elegans software that allows for greater flexibility and accuracy than our prior efforts.
""",
    width = 640,
    justification = :left,
    word_wrap = true,
    fontsize = 10,
)

mipav_screenshot = load("mipav.png")
mipav_axis = Axis(
    left_panel[3,1],
    aspect = DataAspect(),
    yreversed = true,
    height = 360,
    title = "Janelia MIPAV allows for deformable worm cross sections"
)
image!(mipav_axis, mipav_screenshot')
hidedecorations!(mipav_axis)
hidespines!(mipav_axis)


let slider = Slider(left_panel[5,1], range = 0:100)
    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
#    mesh!(ax, _mesh; colorrange, color = _colorpha = 0.1)


    twisted_model_axis = Axis3(
        left_panel[4,1],
        height = 300,
        aspect = :data,
        #viewmode = :stretch,
        title = "Twisted Model Segmented via MIPAV"
    )
    model_mesh = ShroffCelegansModels.get_model_contour_mesh(model; transform_points = scale_dims)
    mesh!(twisted_model_axis, model_mesh, color = _color, alpha = 0.8, shading = MakieCore.automatic)

end

#center_ax = Axis(center_panel[1,1])
center_box = Box(
    center_panel[1:2,1],
    color = colors[2],
    height = 972,
    cornerradius=10
)
center_title = Label(
    center_panel[1,1],
    "Modeling, untwisting and warping the C. elegans embryo in 4D",
    fontsize = 48,
    width = 640,
    word_wrap = true,
    color = :white
    #backgroundcolor = colors[3]
)
center_authors = Label(
    center_panel[2,1],
    fontsize = 24,
    width = 600,
    word_wrap = true,
    color = :white,
    "Mark Kittisopikul, Ryan Christensen, Grant Kroeschell, Johnny Bui, Tosif Ahamed, Hari Shroff"
)


logo_axis = Axis(right_panel[1,1], aspect = DataAspect(), yreversed = true, height = 180)
#logo_axis = Axis3(right_panel[1,1])
scicomp_soft_logo = load("scicomp_soft_logo.png")
scicomp_soft_logo_out = image!(logo_axis, scicomp_soft_logo')
hidedecorations!(logo_axis)
hidespines!(logo_axis)

function poster_animate_untwist(smodel, layout)
    model = smodel.twisted_model
    _pts, _faces = ShroffCelegansModels.get_model_manifold_mesh_components(model)
    r = LinRange(0,1,length(smodel))
    spts = ShroffCelegansModels.central_spline(smodel).(r)

    # fig = Figure()

    # ax = Axis3(fig[1,1], aspect=:data)

    cs_pts = ShroffCelegansModels.central_spline(model).(r)
    scs_pts = ShroffCelegansModels.central_spline(smodel).(r)

    g = GridLayout(layout)

    out = lines(g[1,1], [cs_pts; swapyz.(scs_pts)]; visible = false)
    ax = out.axis

    M, pts = stretch_worm(_pts, _faces, spts, 0.0)
    M = Observable(M)
    L = Observable(@view(pts[1:3:end]))
    C = Observable(@view(pts[2:3:end]))
    R = Observable(@view(pts[3:3:end]))

    model_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points = p->swapyz(p - cs_pts[1])))
    mesh_plot = mesh!(ax, model_mesh; transparency = true)

    mesh!(ax, M; shading = MakieCore.automatic, color = [i for c in 1:length(model) for i in 1:3], colormap = :buda)
    lines!(ax, L, color=:red)
    lines!(ax, C, color=:magenta)
    lines!(ax, R, color=:green)

    s = Slider(g[2,1], range=0:31)

    on(s.value) do value

        if value < 10
            mesh_plot.alpha[] = (100-value*10)/100
        else
            mesh_plot.alpha[] = 0
        end

        if value <= 20
            M[], pts = stretch_worm(_pts, _faces, spts, 0.1*(max(value,10)-10))
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
        end

        if value >= 21 && value <= 31
            M[], pts = untwist_worm(_pts, _faces, spts, 0.1*(value-21))
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
        end
    end
end

Label(
    right_panel[2,1],
    text = "Untwisting Process",
    width = 480,
    font = :bold,
    justification = :center
)

poster_animate_untwist(smodel, right_panel[3,1])


#=
untwisted_model_axis = Axis3(
    right_panel[4,1],
    height = 360,
    aspect = :data,
    #viewmode= :stretch,
    title = "Straightened Model"
)
smodel_mesh = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points = swapyz_scale)
mesh!(untwisted_model_axis, smodel_mesh)
=#

function poster_average_models_with_annotations(
    layout,
    models,
    dataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    #f = Figure(size = (1920, 1080))
    title = Observable("Title")

    f = GridLayout(layout)

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
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom), visible = false)
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
        title[] = "Average t = $value"
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

end

poster_average_models_with_annotations(right_panel[4,1],
avg_models, int_ds[1]; use_myuntwist=true)

poster_figure = f

display(poster_figure)