function meshscatter_all()
    fig = Figure(size = (1920, 1080))
    ax = LScene(fig[1,1])
    coordinates = values(my_annotation_position_cache)
    #T = MeshScatter{Tuple{Vector{Point{3, Float64}}}}
    T = Observable{Vector{Point{3, Float64}}}
    s = Vector{T}(undef, length(coordinates))
    _markersize = Observable(0.1)
    for (i, v) in enumerate(coordinates)
        s[i] = Observable(v[end])
        meshscatter!(s[i], markersize=_markersize)
    end
    display(fig)
    time_points = axes(first(coordinates), 1)
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
    fig
end
#with_theme(meshscatter_all, theme_black())
set_theme!(theme_black())