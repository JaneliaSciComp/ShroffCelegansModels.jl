using Makie
using Makie: throttle, Button
using Printf
using GeometryBasics

function meshscatter_average(average_annotations_dict; nerve_ring = false, models = nothing)
    fig = Figure(size = (1920, 1080))
    ax = LScene(fig[1,1:4]; show_axis = false)
    coordinates = values(average_annotations_dict)
    #T = MeshScatter{Tuple{Vector{Point{3, Float64}}}}
    T = Observable{Vector{Point{3, Float64}}}
    s = Vector{T}(undef, length(coordinates))
    _markersize = Observable(0.5)
    _keys = collect(keys(average_annotations_dict))

    colors_dict = load_colors_dict()
    function get_color(annotation)
        annotation = lowercase(annotation)
        annotation = replace(annotation, "/" => "_")
        get(colors_dict, annotation, RGBf(1,1,1))
    end

    if nerve_ring
        nerve_ring_data = average_annotations_dict["DCR6485_RPM1_NU"]
        p = sortperm(nerve_ring_data.annotations)
        nerve_ring_positions = Observable(nerve_ring_data.positions[end][p])
        nerve_ring_line = lines!(nerve_ring_positions, color = :blue, linewidth = 2)
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v.positions[end])
            if contains(_keys[i], "DCR6485_RPM1_NU")
                colors = get_color.(v.annotations)
                function inspector_label(self, idx, position)
                    return string(v.annotations[idx], position)
                end
                meshscatter!(s[i], markersize=_markersize, color = colors, inspector_label = inspector_label)
            else
                meshscatter!(s[i], markersize=_markersize, color = :grey, alpha = 0.1)
            end
        end
    else
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v.positions[end])
            colors = get_color.(v.annotations)
            function inspector_label(self, idx, position)
                return string(v.annotations[idx], position)
            end
            meshscatter!(s[i], markersize=_markersize, color = colors, inspector_label = inspector_label)
        end
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
    time_text = Observable("hpf = 14:00")
    text!(-20, 0, 20; text = time_text)
    text!( 20, 190, -20; text = "10 μm")
    scalebar = Observable(Point3f[[21, 190, -21], [21, 200, -21]])
    lines!(scalebar, color = :white, linewidth = 5)
    cc = Camera3D(ax.scene;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(60, 90, 0)
    )

    time_points = axes(first(coordinates).positions, 1)
    time_slider = Makie.Slider(fig[2,1:4], range = time_points, startvalue = 201)
    #grid = GridLayout(tellwidth = false, tellheight = false, height = Fixed(5))
    xy_button = Button(fig; label = "XY", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    xz_button = Button(fig; label = "XZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    yz_button = Button(fig; label = "YZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    markersize_menu = Menu(fig, options = string.(0.5:0.1:4), default = "0.5",
        cell_color_inactive_even = RGBf(0.5, 0.5, 0.5),
        cell_color_inactive_odd = RGBf(0.5, 0.5, 0.5),
    )
    #grid[1,1] = button
    fig[3,1] = xy_button
    fig[3,2] = xz_button
    fig[3,3] = yz_button
    fig[3,4] = markersize_menu
    zoom!(ax.scene, 2)
    on(throttle(0.1, time_slider.value)) do t
        total_minutes = (t-1)/200*420
        hours = 7 + round(Int, total_minutes/60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        time_text[] = "hpf = $hours:$(@sprintf("%02d", minutes))"
        for (i, v) in enumerate(coordinates)
            s[i][] = v.positions[t]
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
        scalebar[] = Point3f[[21, 190, -21], [11, 190, -21]]
    end
    on(xz_button.clicks) do _
        update_cam!(ax.scene, cc, 0, π/2)
        scalebar[] = Point3f[[21, 190, -21], [21, 200, -21]]
    end
    on(yz_button.clicks) do _
        update_cam!(ax.scene, cc, 0, 0)
        scalebar[] = Point3f[[21, 190, -21], [21, 200, -21]]
    end
    on(markersize_menu.selection) do _
        _markersize[] = parse(Float64, markersize_menu.selection[])
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
    fig
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
