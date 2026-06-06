#julia> fig = with_theme(theme_black()) do
using Colors
using TiffImages
using FixedPointNumbers
function build_figure(model=models[49], vol = nothing, redch = nothing)
    r = LinRange(0,1,length(model))
    verts = map(model.transverse_splines[[1,9,17,25]]) do s
           s.(r)
    end |> stack
    central_pts = model.central_spline.(r);
    function square(n)
        pts = [
            verts[n, 1] + verts[n, 2],
            verts[n, 3] + verts[n, 2],
            verts[n, 3] + verts[n, 4],
            verts[n, 1] + verts[n, 4]
        ] .- central_pts[n]
        # println(pts)
        GeometryBasics.Mesh(pts, [1, 2, 3, 3, 4, 1])
    end
    smodel = ShroffCelegansModels.Types.StraightenedCelegansModel(model)
    _square_obs = Observable(square(1))
    _contour_obs = Observable(ShroffCelegansModels.get_model_contour_mesh(model))


    green_colormap = to_colormap(:greens)
    green_colormap[1] = RGBA(0,0,0,0)
    red_colormap = to_colormap(:reds)
    red_colormap[1] = RGBA(0,0,0,0)
    blue_colormap = to_colormap(:blues)
    blue_colormap[1] = RGBA(0,0,0,0)
    
    sections = ShroffCelegansModels.get_sections(model)
    straight_sections = ShroffCelegansModels.get_sections(smodel)
    
    _straight_contour_obs = Observable(
        ShroffCelegansModels.get_model_contour_mesh(
            straight_sections,
            transform_points = swapyz_scale,
        )
    )

    #vol = TiffImages.load(raw"X:\shrofflab\RW10598\Data\598_Slitscan_6um_5min_Pos1\Decon_registered\RegB\Decon_reg_67.tif");
    if isnothing(vol)
        vol = TiffImages.load(raw"X:\shrofflab\OD1599_NU\112719_Pos3\Decon_Reg\RegB\Decon_reg_60.tif")
        redch = TiffImages.load(raw"X:\shrofflab\OD1599_NU\112719_Pos3\Decon_Reg\RegA\Decon_reg_60.tif")
    end
    raw_vol = reinterpret(N0f32, vol)

    fig = Figure()
    ax = LScene(fig[1:2, 1]; show_axis=true)
    t = 1:length(model)
    ax_straight = LScene(fig[3,1]; show_axis=true)
    sliders = SliderGrid(fig[4, 1],
        (label="Spline Parameter", range=t),
    )
    r = LinRange(0,1,length(model))
    low_raw_vol = similar(raw_vol)
    low_raw_vol .= 0
    mask = raw_vol .< 0.01
    low_raw_vol[mask] .= raw_vol[mask]
    #v = volume!(ax, permutedims(vol, (2, 1, 3)); colormap=green_colormap, colorrange=(0.00, 0.5))
    #v2 = volume!(ax, permutedims(reinterpret(N0f16, redch), (2, 1, 3)); colormap=red_colormap, colorrange=(0.00, 0.03))
    v = volume!(ax, permutedims(raw_vol, (2,1,3)); colormap=green_colormap, colorrange=extrema(raw_vol))
    #v2 = volume!(ax, permutedims(redch, (2,1,3)); colormap=red_colormap, colorrange=extrema(redch))
    #v_surface = volume!(ax, permutedims(low_raw_vol, (2,1,3)); colormap=blue_colormap, colorrange=(0.00, 0.01), alpha = 0.5)
    mesh!(ax, _contour_obs, alpha=0.5, transparency=true; color = :grey)
    #_square_obs = Observable(square(1))
    mesh!(ax, _square_obs, color=:red)
    right_spline = Observable(ShroffCelegansModels.Types.transverse_spline(model, 1).(r))
    left_spline = Observable(ShroffCelegansModels.Types.transverse_spline(model, 17).(r))
    central_spline = Observable(ShroffCelegansModels.Types.central_spline(model).(r))
    lines!(ax, left_spline, color=:red)
    lines!(ax, right_spline, color=:green)
    lines!(ax, central_spline, color=:magenta)
    pts = seam_cell_pts(model, 0)
    scatter!(ax, pts, color=:red)

    mesh!(ax_straight, _straight_contour_obs, alpha = 0.5, transparency=true; color = :grey)


    #scatter!(ax, pts, color=:red)
    cc = Camera3D(ax.scene)

    zoom!(ax.scene, cc, 0.5)

    cc_straight = Camera3D(ax_straight.scene;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(30, 90, 0)
    )

    zoom!(ax_straight.scene, cc_straight, 0.5)
    update_cam!(ax_straight.scene, cc_straight)

    on(sliders.sliders[1].value) do i
        r = LinRange(0,1,length(model))[1:i]
        _square_obs[] = square(i)
        _contour_obs[] = ShroffCelegansModels.get_model_contour_mesh(
            sections[1:i]; ellipse_points=32)
        right_spline[] = ShroffCelegansModels.Types.transverse_spline(model, 1).(r)
        left_spline[] = ShroffCelegansModels.Types.transverse_spline(model, 17).(r)
        central_spline[] = ShroffCelegansModels.Types.central_spline(model).(r)

        _straight_contour_obs[] = ShroffCelegansModels.get_model_contour_mesh(
            straight_sections[1:i];
            ellipse_points=32,
            transform_points=swapyz_scale
        )
    end

    (;
        fig,
        ax,
        ax_straight,
        cc,
        cc_straight
    )

    Foo(fig) do fig
        display(fig, update=false)
        zoom!(ax.scene, cameracontrols(ax), 0.5)
        zoom!(ax_straight.scene, cameracontrols(ax_straight), 0.2)
    end
end

#=
julia> record(fig, "2024_10_24_sweep_v1.mp4", 1:length(models[60]); visible=true) do i
           _square_obs[] = square(i)
           _contour_obs[] = ShroffCelegansModels.get_model_contour_mesh(sections[1:i]; ellipse_points=32)
           sleep(0.01)
       end
=#

function square_samples(n, radius=nothing)
    if isnothing(radius)
        pts = [
            verts[n, 1] + verts[n, 2],
            verts[n, 3] + verts[n, 2],
            verts[n, 3] + verts[n, 4],
            verts[n, 1] + verts[n, 4]
        ] .- central_pts[n]
    else
        up_vec = normalize(verts[n, 2] .- central_pts[n])
        right_vec = normalize(verts[n, 1] .- central_pts[n])
        pts = [
            up_vec + right_vec,
            up_vec - right_vec,
            -up_vec - right_vec,
            -up_vec + right_vec
        ]
        # pts = normalize.(pts)
        pts .*= radius
        pts .+= (central_pts[n],)
    end
    w = floor(Int, norm(pts[1] - pts[4]))
    h = floor(Int, norm(pts[2] - pts[1]))
    # @info w h norm(up_vec) norm(right_vec) norm(pts[1] - pts[4]) norm(pts[2] - pts[1]) radius
    left = LinRange(pts[4], pts[3], h)
    top = LinRange(Point3f(0), pts[1] - pts[4], w)
    # voxels = Matrix{Point3f}(undef, h, w)
    voxels = Point3f[]
    for start in left
        for pt in (start,) .+ top
            push!(voxels, pt)
        end
    end
    return reshape(voxels, (w, h))
end