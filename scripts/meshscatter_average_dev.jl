using Makie
using Makie: throttle, Button
using Printf
using GeometryBasics

include("crop_video.jl")

function meshscatter_average(average_annotations_dict; nerve_ring = false, models = nothing, session = nothing, xy_bounding_radius = -1)
    max_columns = 1

    fig = Figure(size = (1920, 1080))
    ax = LScene(fig[1,1:max_columns]; show_axis = false)
    coordinates = values(average_annotations_dict)
    #T = MeshScatter{Tuple{Vector{Point{3, Float64}}}}
    T = Observable{Vector{Point{3, Float64}}}
    s = Vector{T}(undef, length(coordinates))
    _markersize = Observable(1.0)
    _keys = collect(keys(average_annotations_dict))

    colors_dict = load_colors_dict()
    function get_color(annotation)
        annotation = lowercase(annotation)
        annotation = replace(annotation, "/" => "_")
        get(colors_dict, annotation, RGBAf(0,0,0,0))
    end

    # Font size
    _fontsize = Observable(40)

    # HPF Label
    label_offset = 8
    time_text = Observable("hpf = 14:00")
    text!(-label_offset, 0, label_offset; text = time_text, fontsize=_fontsize)

    # Scalebar
    scalebar_size_um = 10
    scalebar_y_offset = 180
    scalebar_y_offset_obs = Observable(180)
    scalebar_text = Observable("$scalebar_size_um μm")
    scalebar = Observable(Point3f[
        [label_offset+1, scalebar_y_offset,                    -label_offset-1],
        [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1]
    ])
    scalebar_label_position = lift(scalebar) do pos
        first(pos)
    end
    text!( scalebar_label_position; text = scalebar_text, fontsize = _fontsize, align = (:left, :bottom))
    lines!(scalebar, color = :white, linewidth = 5)
    alpha_obs = Observable(1.0)

    if nerve_ring
        alpha_obs[] = 0.2
        _other_markersize = lift(x->x/2, _markersize)
        nerve_ring_data = average_annotations_dict["DCR6485_RPM1_NU"]
        p = sortperm(nerve_ring_data.annotations)
        nerve_ring_positions = Observable(nerve_ring_data.positions[end][p])
        nerve_ring_line = lines!(nerve_ring_positions, color = :blue, linewidth = 3)
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v.positions[end])
            function inspector_label(self, idx, position)
                return string(v.annotations[idx], position)
            end
            if contains(_keys[i], "DCR6485_RPM1_NU")
                colors = get_color.(v.annotations)
                colors = map(colors) do color
                    return RGBf(0,0,1)
                end
                meshscatter!(s[i], markersize=_markersize, color = colors, inspector_label = inspector_label)
            else
                colors = map(get_color.(v.annotations)) do color
                    if color == RGBAf(0,0,0,0)
                        return color
                    else
                        return RGBf(1, 1, 1)
                    end
                end
                meshscatter!(s[i], markersize=_other_markersize, color = colors, alpha = alpha_obs, inspector_label = inspector_label)
            end
        end
    else
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v.positions[end])
            colors = get_color.(v.annotations)
            function inspector_label(self, idx, position)
                return string(v.annotations[idx], position)
            end
            meshscatter!(s[i], markersize=_markersize, color = colors, inspector_label = inspector_label, alpha = alpha_obs)
        end
        scatter_legend_pts = Point3f[
            [0,  0, 0], 
            [0,  1, 0], 
            [0,  2, 0], 
            [0,  3, 0], 
            [0,  4, 0], 
            [0,  5, 0], 
            [0,  6, 0], 
            [0,  7, 0], 
            [0,  8, 0], 
        ]
        # Horizontal offset between values
        scatter_legend_pts .*= Point3f(0, 20, 0)
        # Vertical offset from worm
        scatter_legend_pts .+= Point3f(label_offset+3, 0, -label_offset-3)
        legend_colors = RGBf[
            RGBf(185,   0, 255),
            RGBf(  0, 255, 255),
            RGBf(255,   0,   0),
            RGBf(255,  99,  93),
            RGBf(  0,   0, 128),
            RGBf(255, 251,  93),
            RGBf( 16, 185,   0),
            RGBf(  0, 128,   0),
            RGBf( 93, 255, 109),
        ]./255
        legend_text = String[
            "Neurons",
            "Neuroblast",
            "Muscle",
            "Intestine",
            "Nerve Ring",
            "Glial",
            "Seam Cell",
            "Other Hypodermal",
            "Pharyngeal",
        ]
        meshscatter!(
            scatter_legend_pts,
            markersize = _markersize,
            color = legend_colors,
            alpha = alpha_obs
        )
        text!(
            scatter_legend_pts;
            text = legend_text,
            align = (:left, :center),
            offset = (10, 0)
        )
    end
    if !isnothing(models)
        seam_cells = Observable(swapyz_scale.(seam_cell_pts(models[end], 2)))
        seam_cell_names = models[end].names
        seam_cell_names[1:2:end] .= replace.(seam_cell_names[1:2:end], 'L' => 'R')
        seam_cell_names = [seam_cell_names[1:2:end]; seam_cell_names[2:2:end]]
        colors = get_color.(seam_cell_names)
        function _inspector_label(self, idx, position)
            return string(seam_cell_names[idx], position)
        end
        meshscatter!(seam_cells, markersize=_markersize, color = colors, inspector_label = _inspector_label)
    end

    cc = Camera3D(ax.scene;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(60, 90, 0)
    )

    time_points = axes(first(coordinates).positions, 1)
    controls_visible = Observable(true)
    time_slider = Makie.Slider(fig[2,1:max_columns], range = time_points, startvalue = 201)
    #grid = GridLayout(tellwidth = false, tellheight = false, height = Fixed(5))
    xy_button = Button(fig; label = "XY", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    xz_button = Button(fig; label = "XZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    yz_button = Button(fig; label = "YZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    markersize_label = Label(fig, text = "Marker Size")
    markersize_menu = Menu(fig, options = string.(0.5:0.1:4), default = "1.0",
        cell_color_inactive_even = RGBf(0.5, 0.5, 0.5),
        cell_color_inactive_odd = RGBf(0.5, 0.5, 0.5),
    )
    record_button = Button(fig; label = "Record", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    record_button2 = Button(fig; label = "Record Loop", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    fontsize_label = Label(fig, text = "Font Size")
    fontsize_menu = Menu(fig, options = string.(20:60), default=string(_fontsize[]))
    transparency_label = Label(fig, text = "Transparency")
    transparency_slider = Makie.Slider(fig, range=0:0.01:1, startvalue = alpha_obs[])
    #grid[1,1] = button
    g = fig[3,1] = GridLayout()
    g[1,1] = xy_button
    g[1,2] = xz_button
    g[1,3] = yz_button
    g[1,4] = markersize_label
    g[1,5] = markersize_menu
    g[1,6] = fontsize_label
    g[1,7] = fontsize_menu
    g[1,8] = transparency_label
    g[1,9] = transparency_slider
    g[1,10] = record_button
    g[1,11] = record_button2
    for c in contents(g)
        c.blockscene.visible = controls_visible
    end
    time_slider.blockscene.visible = controls_visible

    zoom!(ax.scene, 4)
    on(throttle(0.1, time_slider.value)) do t
        total_minutes = (t-1)/200*420
        hours = 7 + round(Int, total_minutes/60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        time_text[] = "hpf = $hours:$(@sprintf("%02d", minutes))"
        if xy_bounding_radius > 0
            for (i, v) in enumerate(coordinates)
                p = v.positions[t]
                map(p) do pt
                    if pt[1].^2 + pt[3].^2 > xy_bounding_radius^2
                        return Point3(NaN)
                    else
                        return pt
                    end
                end
                s[i][] = p
            end
        else
            for (i, v) in enumerate(coordinates)
                s[i][] = v.positions[t]
            end
        end
        if nerve_ring
            nerve_ring_positions[] = nerve_ring_data.positions[t][p]
        end
        if !isnothing(models)
            seam_cells[] = swapyz_scale.(seam_cell_pts(models[t], 2))
        end
    end
    on(xy_button.clicks) do _
        update_cam!(ax.scene, cc, π/2, 0)
        scalebar_y_offset_obs[] = 0
        scalebar[] = Point3f[[label_offset+1, 0, -label_offset-1], [label_offset+1 - scalebar_size_um/2, 0, -label_offset-1]]
        scalebar_text[] = "$(scalebar_size_um÷2) μm"
    end
    on(xz_button.clicks) do _
        update_cam!(ax.scene, cc, 0, π/2)
        scalebar[] = Point3f[[label_offset+1, scalebar_y_offset, -label_offset-1], [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1]]
    end
    on(yz_button.clicks) do _
        update_cam!(ax.scene, cc, 0, 0)
        scalebar[] = Point3f[[label_offset+1, scalebar_y_offset, -label_offset-1], [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1]]
    end
    vid = Observable(DOM.div("Press record..."; id="video_recording", style="color: white; display: none;"))
    on(record_button.clicks) do _
        # record(fig, "/var/www/shroff/test.mp4", time_slider.range[]; update=false) do t
        # VideoStream
        controls_visible[] = false
        vs = Record(fig, time_slider.range[]; update=false) do t
            set_close_to!(time_slider, t)
        end
        vid[] = DOM.div(crop_video(vs); id = "video_recording")
        controls_visible[] = true
        if !isnothing(session)
            evaljs(session, js"""document.getElementById("video_recording").scrollIntoView(true)""")
        end
    end
    on(record_button2.clicks) do _
        # record(fig, "/var/www/shroff/test.mp4", time_slider.range[]; update=false) do t
        # VideoStream
        controls_visible[] = false
        L = time_slider.range[].stop
        _range = 1:2L
        vs = Record(fig, _range; update=false) do t
            if t <= L
                set_close_to!(time_slider, t)
            elseif t <= 2L
                set_close_to!(time_slider, 2L - t)
            end
        end
        vid[] = DOM.div(crop_video(vs); id = "video_recording")
        controls_visible[] = true
        if !isnothing(session)
            evaljs(session, js"""document.getElementById("video_recording").scrollIntoView(true)""")
        end
    end

    on(markersize_menu.selection) do _
        _markersize[] = parse(Float64, markersize_menu.selection[])
    end
    on(fontsize_menu.selection) do _
        _fontsize[] = parse(Int, fontsize_menu.selection[])
    end
    on(transparency_slider.value) do v
        alpha_obs[] = v
    end
    #=
    println("Press any key to continue:")
    readline()
    for i in 1:10
        _markersize[] = 0.1*i
        for t in time_points
            for (i, v) in enumerate(coordinates)
                s[i][] = v[t]
            end
            sleep(0.05)
        end
    end
    =#
    #DOM.body(fig, style=Styles(CSS("background-color" => "black")))
    DataInspector(fig; backgroundcolor = :black)
    DOM.div(fig, vid)
end
#with_theme(meshscatter_all, theme_black())
#set_theme!(theme_black())

function load_colors_dict()
    colors_dict = Dict{String, RGBf}()
    colors_df = CSV.read("colors.csv", DataFrame)
    for row in eachrow(colors_df)
        colors_dict[lowercase(row.Annotation)] = RGBf(row.R/255, row.G/255, row.B/255)
    end
    return colors_dict
end
