using CairoMakie
using ShroffCelegansModels
using ShroffCelegansModels: get_pretwitch_df, get_pretwitch_points_at_time,
    pretwitch_reference_axes, pretwitch_cells_by_name, Point3f, Vec3f, normalize, dot, cross

# Same reconstruction as movie_pretwitch_axes_regression.jl (see that file for
# the full explanation), but with the color scale renormalized EVERY FRAME to
# that frame's own active-cell weight range, instead of held fixed to the
# global max across all 361 frames. Trade-off: more visual contrast within
# any single frame (useful since most individual cells have small weight
# relative to the handful of dominant early lineage nodes -- see the fixed-
# scale version, where most points look pale after those nodes divide away),
# at the cost of the colorbar no longer meaning the same absolute weight
# value from frame to frame.

pretwitch_df = get_pretwitch_df()
ref = pretwitch_reference_axes(pretwitch_df)
grouped = pretwitch_cells_by_name(pretwitch_df)
times = sort(unique(pretwitch_df.time))

println("Precomputing per-frame centroid...")
centroid_of = Dict{Int,Point3f}()
for t in times
    pts = collect(values(get_pretwitch_points_at_time(pretwitch_df, t)))
    centroid_of[t] = Point3f(sum(pts) / length(pts))
end

const MIN_DIST = 10f0
norm3(v) = sqrt(sum(abs2, v))

println("Building raw (un-spliced) direction tracks...")
dir_tracks = Dict{String,Dict{Int,Vector{Float64}}}()
for (name, rows) in grouped
    d = Dict{Int,Vector{Float64}}()
    for (t, p) in rows
        v = p - centroid_of[t]
        norm3(v) < MIN_DIST && continue
        d[t] = Float64.(normalize(v))
    end
    isempty(d) || (dir_tracks[name] = d)
end

println("Loading regression weights...")
csv_path = joinpath(@__DIR__, "..", "pretwitch_axis_regression_weights.csv")
dv_w = Dict{String,Float64}()
lr_w = Dict{String,Float64}()
open(csv_path) do io
    readline(io)  # header
    for line in eachline(io)
        parts = split(line, ',')
        name = parts[1]
        dv_w[name] = parse(Float64, parts[3])
        lr_w[name] = parse(Float64, parts[4])
    end
end
println("$(length(dv_w)) weighted cells loaded.")

const DIVERGING_CMAP = :RdBu
const MIN_RANGE = 1e-3  # floor so a frame with near-zero weights doesn't blow up the color scale

function reconstruct(weights, t)
    pred = zeros(3)
    for (name, w) in weights
        d = get(dir_tracks, name, nothing)
        isnothing(d) && continue
        p = get(d, t, nothing)
        isnothing(p) && continue
        pred .+= w .* p
    end
    nrm = norm3(pred)
    return nrm < 1e-8 ? nothing : Vec3f(pred ./ nrm)
end

frame_data = map(times) do t
    pts_dict = get_pretwitch_points_at_time(pretwitch_df, t)
    names = collect(keys(pts_dict))
    pts = collect(values(pts_dict))
    centroid = centroid_of[t]

    ap = ref.ap
    dv_raw = reconstruct(dv_w, t)
    lr_raw = reconstruct(lr_w, t)
    dv = normalize(dv_raw - dot(dv_raw, ap) * ap)
    lr_candidate = normalize(cross(ap, dv))
    lr = dot(lr_candidate, lr_raw) < 0 ? -lr_candidate : lr_candidate

    stabilized = [Point3f(dot(p - centroid, ap), dot(p - centroid, lr), dot(p - centroid, dv)) for p in pts]
    dv_colors = [get(dv_w, name, 0.0) for name in names]
    lr_colors = [get(lr_w, name, 0.0) for name in names]
    dv_range = max(maximum(abs, dv_colors), MIN_RANGE)
    lr_range = max(maximum(abs, lr_colors), MIN_RANGE)

    (t=t, pts=pts, centroid=centroid, ap=ap, dv=dv, lr=lr, stabilized=stabilized,
     dv_colors=dv_colors, lr_colors=lr_colors, dv_range=dv_range, lr_range=lr_range)
end

fig = Figure(size=(1500, 1300))

ax_raw_dv = Axis3(fig[1, 1], title="raw, colored by dv_weight (per-frame scale)", aspect=:data,
    limits=(50, 350, 20, 220, 40, 240))
ax_stab_dv = Axis3(fig[1, 2], title="stabilized, colored by dv_weight (per-frame scale)", aspect=:data,
    limits=(-150, 150, -150, 150, -150, 150), xlabel="AP", ylabel="LR", zlabel="DV")
ax_raw_lr = Axis3(fig[2, 1], title="raw, colored by lr_weight (per-frame scale)", aspect=:data,
    limits=(50, 350, 20, 220, 40, 240))
ax_stab_lr = Axis3(fig[2, 2], title="stabilized, colored by lr_weight (per-frame scale)", aspect=:data,
    limits=(-150, 150, -150, 150, -150, 150), xlabel="AP", ylabel="LR", zlabel="DV")

pts_obs = Observable(frame_data[1].pts)
stab_obs = Observable(frame_data[1].stabilized)
dv_color_obs = Observable(frame_data[1].dv_colors)
lr_color_obs = Observable(frame_data[1].lr_colors)
dv_colorrange_obs = Observable((-frame_data[1].dv_range, frame_data[1].dv_range))
lr_colorrange_obs = Observable((-frame_data[1].lr_range, frame_data[1].lr_range))

scatter!(ax_raw_dv, pts_obs, color=dv_color_obs, colormap=DIVERGING_CMAP, colorrange=dv_colorrange_obs, markersize=8)
scatter!(ax_stab_dv, stab_obs, color=dv_color_obs, colormap=DIVERGING_CMAP, colorrange=dv_colorrange_obs, markersize=8)
scatter!(ax_raw_lr, pts_obs, color=lr_color_obs, colormap=DIVERGING_CMAP, colorrange=lr_colorrange_obs, markersize=8)
scatter!(ax_stab_lr, stab_obs, color=lr_color_obs, colormap=DIVERGING_CMAP, colorrange=lr_colorrange_obs, markersize=8)

Colorbar(fig[1, 3], colormap=DIVERGING_CMAP, colorrange=dv_colorrange_obs, label="dv_weight (this frame)")
Colorbar(fig[2, 3], colormap=DIVERGING_CMAP, colorrange=lr_colorrange_obs, label="lr_weight (this frame)")

arrow_scale = 60f0
ap_obs = Observable([frame_data[1].centroid, frame_data[1].centroid + arrow_scale * frame_data[1].ap])
dv_obs = Observable([frame_data[1].centroid, frame_data[1].centroid + arrow_scale * frame_data[1].dv])
lr_obs = Observable([frame_data[1].centroid, frame_data[1].centroid + arrow_scale * frame_data[1].lr])

for ax in (ax_raw_dv, ax_raw_lr)
    lines!(ax, ap_obs, color=:black, linewidth=4)
    lines!(ax, dv_obs, color=:black, linewidth=4, linestyle=:dash)
    lines!(ax, lr_obs, color=:black, linewidth=4, linestyle=:dot)
end
for ax in (ax_stab_dv, ax_stab_lr)
    lines!(ax, [Point3f(0,0,0), Point3f(arrow_scale,0,0)], color=:black, linewidth=4)
    lines!(ax, [Point3f(0,0,0), Point3f(0,arrow_scale,0)], color=:black, linewidth=4, linestyle=:dot)
    lines!(ax, [Point3f(0,0,0), Point3f(0,0,arrow_scale)], color=:black, linewidth=4, linestyle=:dash)
end

Legend(fig[1:2, 4],
    [LineElement(color=:black, linewidth=4), LineElement(color=:black, linewidth=4, linestyle=:dash), LineElement(color=:black, linewidth=4, linestyle=:dot)],
    ["AP (Approach 1, fixed)", "DV (regression weighted sum)", "LR (regression weighted sum)"],
    "Legend")

title_obs = Observable("t=$(frame_data[1].t)  (regression-weighted axes, per-frame color scale)")
Label(fig[0, 1:4], title_obs, fontsize=20, tellwidth=false)

snapshot_indices = Set([1, 91, 181, 271, 361])

outpath = joinpath(@__DIR__, "..", "pretwitch_axes_regression_perframe.mp4")
record(fig, outpath, eachindex(frame_data); framerate=24) do i
    d = frame_data[i]
    pts_obs[] = d.pts
    stab_obs[] = d.stabilized
    dv_color_obs[] = d.dv_colors
    lr_color_obs[] = d.lr_colors
    dv_colorrange_obs[] = (-d.dv_range, d.dv_range)
    lr_colorrange_obs[] = (-d.lr_range, d.lr_range)
    ap_obs[] = [d.centroid, d.centroid + arrow_scale * d.ap]
    dv_obs[] = [d.centroid, d.centroid + arrow_scale * d.dv]
    lr_obs[] = [d.centroid, d.centroid + arrow_scale * d.lr]
    title_obs[] = "t=$(d.t)  (regression-weighted axes, per-frame color scale)"
    if i in snapshot_indices
        save(joinpath(@__DIR__, "..", "pretwitch_axes_regression_perframe_frame$(lpad(i, 3, '0')).png"), fig)
    end
end
println("wrote $outpath")
