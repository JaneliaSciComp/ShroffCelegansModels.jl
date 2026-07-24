using CairoMakie
using ShroffCelegansModels
using ShroffCelegansModels: get_pretwitch_df, get_pretwitch_points_at_time, build_pretwitch_central_spline_series,
    pretwitch_reference_axes, Point3f, Vec3f, dot, cross, normalize, norm,
    left_seam_cells, right_seam_cells

# Verification movie for the pretwitch central spline (see
# src/demo_averaging/pretwitch_straighten.jl). The fitted central spline is
# now drawn ONLY through TERMINAL seam cells -- pairs that have reached their
# exact final differentiated identity, not merely a still-dividing shared
# ancestor (see PretwitchPairTrack.is_terminal). No pair is terminal before
# frame 170 (~47% into the pretwitch window), so the curve/knots are simply
# absent for those early frames; the surrounding ancestor-cell dots/labels
# and L-R connectors are still shown throughout as context.
#
# Each dot is one labeled point per DISTINCT currently-active ancestor cell
# (not one per seam cell -- many of the 20 L/R seam-cell ancestors share the
# same actual cell before they've diverged, so plotting all 20 labels would
# stack duplicates on top of each other). Each label is "<actual cell name>
# (<seam cells descending from it>)", e.g. "ABar (H2L,V1L,V2L,V4L,V6L,H2R,
# V1R,V2R,V4R,V6R)". A line is still drawn between each of the 10 canonical
# L/R pairs' current positions (zero-length if that pair hasn't diverged
# yet), alongside the existing Approach-1 (WormGuides compass) AP axis arrow
# for a rough visual reference. The title shows the running terminal-station
# count and whether THIS frame's nominal knot order also happens to be
# spatially monotonic along that AP axis -- a diagnostic only, not something
# that alters the fit.
#
# Renders two outputs via render_pretwitch_central_spline_movie (single
# function, `show_surface` keyword toggles the body-surface tube mesh; the
# expensive per-frame precomputation, including tube meshes, is shared
# between both calls): pretwitch_central_spline.mp4 (curve-only, no surface)
# and pretwitch_central_spline_with_surface.mp4 (includes the body surface).

pretwitch_df = get_pretwitch_df()
series = build_pretwitch_central_spline_series(pretwitch_df)
ref = pretwitch_reference_axes(pretwitch_df)

# Rainbow palette, restored from the earlier version of this movie -- one
# FIXED color per canonical station (H0=1st,...,T=10th), stable across every
# frame (unlike the old "color by current spline-knot group" scheme, which
# reassigned colors as stations merged/split from frame to frame).
const STATION_COLORS = [:red, :orange, :gold, :green, :teal, :blue, :purple, :magenta, :brown, :black]

function sample_spline(frame; n=100)
    isnothing(frame.central_spline) && return Point3f[]
    return [Point3f(frame.central_spline(x)) for x in range(0, 1, length=n)]
end

"""
    finite_diff_tangents(curve)

Unit tangent at each point of `curve` via central differences (forward/
backward at the endpoints) -- a simple numerical stand-in for an analytic
derivative, adequate for orienting cross-section circles on the already-
blended `interp_curve` (which has no single analytic spline of its own).
"""
function finite_diff_tangents(curve::Vector{Point3f})
    n = length(curve)
    map(1:n) do j
        d = j == 1 ? curve[2] - curve[1] :
            j == n ? curve[n] - curve[n-1] :
            curve[j+1] - curve[j-1]
        normalize(Vec3f(d))
    end
end

"""
    tube_mesh(interp_curve, right_vecs)

Build a `GeometryBasics.Mesh` tube around `interp_curve`: a circular cross
section at each sample point, oriented by the local tangent (finite
differences, see [`finite_diff_tangents`](@ref)) and the local `right_vecs[j]`
(already scaled to the desired radius -- see how `right_vec` is blended in
the main per-frame loop below), exactly mirroring `build_celegans_model`'s
posttwitch convention: `normal = normalize(cross(tangent, right)) * radius`,
then `get_circle_points(right, normal, center)` for 32 points/ring, stitched
into quads by the same `get_model_contour_mesh` posttwitch uses.
"""
function tube_mesh(interp_curve::Vector{Point3f}, right_vecs::Vector{Vec3f})
    tangents = finite_diff_tangents(interp_curve)
    sections = map(interp_curve, tangents, right_vecs) do center, tangent, rv
        normal = normalize(cross(tangent, rv)) * norm(rv)
        Point3f.(ShroffCelegansModels.get_circle_points(rv, normal, center))
    end
    return ShroffCelegansModels.get_model_contour_mesh(sections)
end

"""
    unique_ancestor_labels(tracks, i)

One `(point, label, color)` per DISTINCT currently-active ancestor cell among
all 20 L/R seam-cell-ancestor positions at pair-track index `i` -- i.e.
deduplicated by actual cell identity, not by seam-cell name, so stations that
haven't yet diverged (sharing one physical cell) get exactly one label
between them. `label` is `"<cell name> (<comma-separated descendant seam
cells>)"`; `color` is `STATION_COLORS` for whichever canonical station (in
nominal H0..T order) first contributes that cell.
"""
function unique_ancestor_labels(tracks, i)
    seen = Dict{String,Tuple{Point3f,Vector{String},Any}}()
    order = String[]
    for (k, tr) in enumerate(tracks)
        for (name, pt, seam_label) in ((tr.left_name[i], tr.left[i], left_seam_cells[k]),
                                        (tr.right_name[i], tr.right[i], right_seam_cells[k]))
            if haskey(seen, name)
                push!(seen[name][2], seam_label)
            else
                seen[name] = (pt, [seam_label], STATION_COLORS[mod1(k, length(STATION_COLORS))])
                push!(order, name)
            end
        end
    end
    pts = [seen[name][1] for name in order]
    texts = ["$(name) ($(join(seen[name][2], ",")))" for name in order]
    colors = [seen[name][3] for name in order]
    return pts, texts, colors
end

# Raw (unmerged) per-station midpoints -- recompute directly since
# PretwitchCentralSplineFrame only stores the terminal-only knot positions.
tracks = ShroffCelegansModels.pretwitch_pair_tracks(pretwitch_df)

# Interpolated curve: medial axis (straight line along the fixed AP
# direction) at/before the first frame a central spline can be fit, the
# fitted spline itself by the final frame, and a plain linear blend of the
# two (sample-by-sample, at matching arc-length-like parameter positions) in
# between -- i.e. w=0 is pure medial axis, w=1 is pure spline.
const FIRST_SPLINE_FRAME = findfirst(f -> !isnothing(f.central_spline), series.frames)
const T_START = series.frames[FIRST_SPLINE_FRAME].time
const T_END = series.frames[end].time

# Fixed per-station line color, one pair (2 vertices) per canonical station,
# always in nominal H0..T order -- never changes frame to frame, so this can
# be a plain (non-Observable) array.
lr_vertex_colors = reduce(vcat, [[c, c] for c in STATION_COLORS[mod1.(1:length(tracks), length(STATION_COLORS))]])

# Precompute all per-frame plot data.
frame_data = map(1:length(series.frames)) do i
    frame = series.frames[i]
    left_pts = [tr.left[i] for tr in tracks]
    right_pts = [tr.right[i] for tr in tracks]
    # Interleaved (left,right,left,right,...) so linesegments! draws one
    # segment per consecutive pair -- exactly the 10 L/R connecting lines.
    lr_segments = reduce(vcat, [[l, r] for (l, r) in zip(left_pts, right_pts)])
    label_pts, label_texts, label_colors = unique_ancestor_labels(tracks, i)
    all_pts = collect(values(get_pretwitch_points_at_time(pretwitch_df, frame.time)))
    # Medial axis: the fixed Approach-1 AP axis, extended to span the full
    # anterior-posterior extent of the ENTIRE point cloud this frame (not a
    # fixed arbitrary length) -- project every cell onto the axis and take
    # the extremes, anchored at the axis's own fixed origin so the line
    # doesn't jitter frame to frame independent of the point cloud's spread.
    projections = [dot(p - ref.origin, ref.ap) for p in all_pts]
    lo, hi = extrema(projections)
    ap_line = [Point3f(ref.origin + lo * ref.ap), Point3f(ref.origin + hi * ref.ap)]

    curve = sample_spline(frame)
    # Medial axis sampled at the same number of points / parameter positions
    # as the spline, so the two are pointwise-comparable for interpolation.
    medial_samples = [Point3f(ref.origin + (lo + (hi - lo) * x) * ref.ap) for x in range(0, 1, length=100)]
    w = clamp((frame.time - T_START) / (T_END - T_START), 0.0, 1.0)
    interp_curve = isempty(curve) ? medial_samples :
                   [Point3f((1 - w) * m + w * s) for (m, s) in zip(medial_samples, curve)]

    # Cross-section radius/orientation: blend, exactly like the curve itself,
    # between an EARLY vector (half the whole point cloud's left-right extent,
    # pointing in the fixed reference LR direction -- a constant-radius,
    # non-twisting cylinder around the straight medial axis) and a LATE
    # vector (the actual local left-right half-distance/direction at each
    # terminal seam-cell pair, from right_vector_spline). Blending the full
    # vector (not radius and direction separately) keeps this identical in
    # spirit to how interp_curve itself blends full positions.
    lr_projections = [dot(p - ref.origin, ref.lr) for p in all_pts]
    lr_lo, lr_hi = extrema(lr_projections)
    early_right_vec = Vec3f(((lr_hi - lr_lo) / 2) * ref.lr)
    right_vecs = if isnothing(frame.right_vector_spline)
        fill(early_right_vec, 100)
    else
        [Vec3f((1 - w) * early_right_vec + w * Vec3f(frame.right_vector_spline(x))) for x in range(0, 1, length=100)]
    end
    surface = tube_mesh(interp_curve, right_vecs)

    (t=frame.time, all_pts=all_pts, lr_segments=lr_segments, label_pts=label_pts, label_texts=label_texts,
     label_colors=label_colors, knots=frame.knot_positions,
     curve=curve, n_terminal=length(frame.terminal_stations),
     ap_line=ap_line, interp_curve=interp_curve, surface=surface, is_ap_monotonic=frame.is_ap_monotonic)
end

title_str(d) = "t=$(d.t)  terminal stations: $(d.n_terminal)/10  spatially AP-monotonic: $(d.is_ap_monotonic)"

"""
    render_pretwitch_central_spline_movie(frame_data; show_surface=true, outpath, snapshot_indices=Set{Int}())

Render the pretwitch central-spline verification movie from precomputed
`frame_data`. `show_surface` toggles the body-surface tube mesh (see
[`tube_mesh`](@ref)) on or off -- everything else (background cells, L/R
connectors, ancestor labels, fitted spline, medial axis, interpolated curve)
is always shown. Set `show_surface=false` for a lighter-weight, curve-only
view when the surface itself isn't the point.
"""
function render_pretwitch_central_spline_movie(frame_data; show_surface::Bool=true, outpath, snapshot_indices=Set{Int}())
    fig = Figure(size=(1300, 800))
    ax = Axis3(fig[1, 1], title="pretwitch central spline", aspect=:data,
        limits=(50, 350, 20, 220, 40, 240))

    all_pts_obs = Observable(frame_data[1].all_pts)
    lr_segments_obs = Observable(frame_data[1].lr_segments)
    label_pts_obs = Observable(frame_data[1].label_pts)
    label_texts_obs = Observable(frame_data[1].label_texts)
    label_colors_obs = Observable(frame_data[1].label_colors)
    knots_obs = Observable(frame_data[1].knots)
    curve_obs = Observable(frame_data[1].curve)
    ap_line_obs = Observable(frame_data[1].ap_line)
    interp_curve_obs = Observable(frame_data[1].interp_curve)

    # Background layer (drawn first, so the highlighted seam-cell points/labels
    # render on top of it): every other pretwitch cell at this frame, for context.
    scatter!(ax, all_pts_obs, color=(:gray, 0.35), markersize=5)
    surface_obs = nothing
    if show_surface
        surface_obs = Observable(frame_data[1].surface)
        mesh!(ax, surface_obs, color=(:skyblue, 0.25), transparency=true)
    end
    linesegments!(ax, lr_segments_obs, color=lr_vertex_colors, linewidth=2)
    scatter!(ax, label_pts_obs, color=label_colors_obs, marker=:circle, markersize=16, strokecolor=:black, strokewidth=1)
    text!(ax, label_pts_obs; text=label_texts_obs, fontsize=11, align=(:left, :bottom), offset=(6, 6))
    scatter!(ax, knots_obs, color=:black, marker=:star5, markersize=22, strokecolor=:white, strokewidth=1)
    lines!(ax, curve_obs, color=:black, linewidth=3)
    lines!(ax, ap_line_obs, color=:red, linewidth=3, linestyle=:dash)
    lines!(ax, interp_curve_obs, color=:dodgerblue, linewidth=4)

    legend_elements = Any[MarkerElement(color=(:gray, 0.35), marker=:circle, markersize=10)]
    legend_labels = String["other pretwitch cells (context)"]
    if show_surface
        push!(legend_elements, PolyElement(color=(:skyblue, 0.25)))
        push!(legend_labels, "body surface (circular cross-sections around interpolated curve)")
    end
    append!(legend_elements, [
        MarkerElement(color=:gray, marker=:circle, markersize=14),
        LineElement(color=:gray, linewidth=2), MarkerElement(color=:black, marker=:star5, markersize=18),
        LineElement(color=:black, linewidth=3), LineElement(color=:red, linewidth=3, linestyle=:dash),
        LineElement(color=:dodgerblue, linewidth=4)])
    append!(legend_labels, [
        "distinct ancestor cell: \"name (descendant seam cells)\" (color = station of first descendant, H0..T)",
        "L–R pair connector (color = station, H0..T)", "spline knot (terminal seam cell, fixed nominal order)",
        "fitted central spline (terminal seam cells only)", "medial axis (Approach 1 AP direction, spanning full point cloud)",
        "interpolated curve (medial axis at t=$(T_START) -> spline at t=$(T_END), linear blend)"])
    Legend(fig[1, 2], legend_elements, legend_labels, "Legend")

    title_obs = Observable(title_str(frame_data[1]))
    Label(fig[0, 1:2], title_obs, fontsize=20, tellwidth=false)

    snapshot_prefix = splitext(outpath)[1]
    record(fig, outpath, eachindex(frame_data); framerate=24) do i
        d = frame_data[i]
        all_pts_obs[] = d.all_pts
        lr_segments_obs[] = d.lr_segments
        label_pts_obs[] = d.label_pts
        label_texts_obs[] = d.label_texts
        label_colors_obs[] = d.label_colors
        knots_obs[] = d.knots
        curve_obs[] = d.curve
        ap_line_obs[] = d.ap_line
        interp_curve_obs[] = d.interp_curve
        show_surface && (surface_obs[] = d.surface)
        title_obs[] = title_str(d)
        if i in snapshot_indices
            save("$(snapshot_prefix)_frame$(lpad(i, 3, '0')).png", fig)
        end
    end
    println("wrote $outpath")
end

snapshot_indices = Set([1, 31, 91, 151, 181, 271, 361])

render_pretwitch_central_spline_movie(frame_data; show_surface=false,
    outpath=joinpath(@__DIR__, "..", "pretwitch_central_spline.mp4"))
render_pretwitch_central_spline_movie(frame_data; show_surface=true,
    outpath=joinpath(@__DIR__, "..", "pretwitch_central_spline_with_surface.mp4"), snapshot_indices)
