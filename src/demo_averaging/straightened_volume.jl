using TiffImages
using FixedPointNumbers
using ShroffCelegansModels
using Interpolations
using Colors

function get_regb(dataset, idx)
    tp = range(dataset.cell_key)[idx]
    joinpath(dataset.path, "Decon_reg_$(tp).tif")
end
idx = 80

dataset = datasets["OD1599_NU"][2]
#vol = TiffImages.load(raw"X:\shrofflab\OD1599_NU\112719_Pos3\Decon_Reg\RegB\Decon_reg_60.tif")
vol = TiffImages.load(get_regb(dataset, idx))
raw_vol = reinterpret(N0f32, vol);
permuted_vol = permutedims(raw_vol, (2, 1, 3))
itp = Interpolations.interpolate(axes(permuted_vol), permuted_vol, Gridded(Interpolations.Linear()))
models, smodels = build_models_over_time(dataset);
r = LinRange(0, 1, length(models[idx]))
max_radius = ShroffCelegansModels.max_radius_function(models[idx]).(r) |> maximum
verts = map(models[idx].transverse_splines[[1, 9, 17, 25]]) do s
    s.(r)
end |> stack;
central_pts = models[idx].central_spline.(r)

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
        ] .* radius .+ (central_pts[n],)
    end
    w = floor(Int, norm(pts[1] - pts[4]))
    h = floor(Int, norm(pts[2] - pts[1]))
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

L = length(models[idx])
pts = square_samples.(1:L, (max_radius,)) |> stack

f(pt) =
    try
        itp(pt...)
    catch err
        NaN
    end

obs_vol = Observable(PermutedDimsArray(f.(@view pts[1:8:end, 1:8:end, 1:8:end]), (2, 3, 1)))
green_colormap = to_colormap(:greens)
green_colormap[1] = RGBA(0,0,0,0)

with_theme(theme_black()) do
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=false)
    volume!(ax, (-48, 48), (1, L), (-48, 48), obs_vol, colorrange=(0.23, 0.27), colormap=green_colormap)
    # volume!(ax, (-48,48),(1,L),(-48,48), obs_volA, colorrange=(0.23,1.0), colormap=red_colormap) 
    cc = Camera3D(ax.scene, projection=:orthographic)
    display(fig, update=true)
    obs_alpha = Observable(0.0)
    mesh!(ShroffCelegansModels.get_model_contour_mesh(smodels[idx], transform_points=swapyz), alpha=obs_alpha, color=:grey)
    record(fig, "rotate_straightened2.mp4", LinRange(-5π / 2, π / 2, 60)) do θ
        update_cam!(ax.scene, cc, 0, θ)
        obs_alpha[] = (θ + π / 2) / π * 0.5
        display(fig, update=false)
        sleep(0.01)
    end
end