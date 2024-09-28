using Makie
using Makie: throttle, Button
using Printf
using GeometryBasics

function meshscatter_all(; nerve_cord = false)
    fig = Figure(size = (1920, 1080))
    ax = LScene(fig[1,1:3]; show_axis = false)
    coordinates = values(my_annotation_position_cache)
    #T = MeshScatter{Tuple{Vector{Point{3, Float64}}}}
    T = Observable{Vector{Point{3, Float64}}}
    s = Vector{T}(undef, length(coordinates))
    _markersize = Observable(0.5)
    _keys = collect(keys(my_annotation_position_cache))
    if nerve_cord
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v[end])
            if contains(_keys[i], "DCR6485_RPM1_NU")
                meshscatter!(s[i], markersize=_markersize, color = :blue)
            else
                meshscatter!(s[i], markersize=_markersize, color = :grey, alpha = 0.1)
            end
        end
    else
        for (i, v) in enumerate(coordinates)
            s[i] = Observable(v[end])
            meshscatter!(s[i], markersize=_markersize)
        end
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

    time_points = axes(first(coordinates), 1)
    time_slider = Makie.Slider(fig[2,1:3], range = time_points, startvalue = 201)
    #grid = GridLayout(tellwidth = false, tellheight = false, height = Fixed(5))
    xy_button = Button(fig; label = "XY", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    xz_button = Button(fig; label = "XZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    yz_button = Button(fig; label = "YZ", buttoncolor = RGBf(0.5, 0.5, 0.5), tellwidth = false)
    #grid[1,1] = button
    fig[3,1] = xy_button
    fig[3,2] = xz_button
    fig[3,3] = yz_button
    zoom!(ax.scene, 2)
    on(throttle(0.1, time_slider.value)) do t
        total_minutes = (t-1)/200*420
        hours = 7 + round(Int, total_minutes/60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        time_text[] = "hpf = $hours:$(@sprintf("%02d", minutes))"
        for (i, v) in enumerate(coordinates)
            s[i][] = v[t]
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
    DOM.body(fig, style=Styles(CSS("background-color" => "black")))
end
#with_theme(meshscatter_all, theme_black())
#set_theme!(theme_black())
