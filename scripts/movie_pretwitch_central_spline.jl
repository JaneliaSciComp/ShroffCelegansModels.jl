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
    points_near_ap_plane(pts, ref, ap_position; min_points=2, initial_thresh=40.0, growth=1.5, max_thresh=300.0)

Cells within `initial_thresh` of the cross-sectional plane through
`ref.origin + ap_position * ref.ap`, perpendicular to `ref.ap` (i.e.
`|dot(p - ref.origin, ref.ap) - ap_position| < thresh`) -- a LOCAL slice at
a specific point along the AP axis, not the whole point cloud. `ap_position`
is a signed scalar offset from `ref.origin` along `ref.ap` (the same units
as `dot(p - ref.origin, ref.ap)`), so evaluating this at many positions along
the AP extent gives a genuinely varying (tapering near the ends of an ovoid
embryo) cross-sectional profile rather than one constant body-wide radius.

`initial_thresh=40` (widened from an initial `15`, which gave a visibly
jagged radius profile -- adjacent 100-sample slices barely overlapped, so a
single point entering/leaving a narrow window caused a sharp jump; measured
directly: max sample-to-sample jump ~16-26 units at thresh=15 vs ~7-16 at
thresh=40-45 across several frames). Wider slices trade a softer, less sharp
taper right at the very tips for much smoother variation in between --
consistent with a real ovoid body rather than sample noise. Callers should
still pass a per-frame-adjusted `initial_thresh` (see
[`adaptive_ap_plane_thresh`](@ref)) rather than this flat default, since a
FIXED width is itself too narrow relative to how sparse early frames are
(as few as 4 cells total at t=0) -- the same 40 units that works well once
hundreds of cells exist is comparatively much sparser-sampled early on.

Very early frames have only a handful of cells total (and slices near the
AP extremes are inherently sparse even later), so the threshold grows (up to
`max_thresh`) until at least `min_points` fall in the slice -- `min_points`
defaults to just 2 (the minimum for a non-degenerate extent) rather than
requiring a fuller sample, since demanding many points near a tapering tip
would defeat the point of a local measurement.
"""
function points_near_ap_plane(pts, ref, ap_position::Real; min_points=2, initial_thresh=40.0, growth=1.5, max_thresh=300.0)
    thresh = initial_thresh
    local slice
    while true
        slice = filter(p -> abs(dot(p - ref.origin, ref.ap) - ap_position) < thresh, pts)
        (length(slice) >= min_points || thresh > max_thresh) && break
        thresh *= growth
    end
    return slice
end

"""
    adaptive_ap_plane_thresh(n_cells, n_reference; base_thresh=40.0, power=0.2)

Per-FRAME slice threshold for [`points_near_ap_plane`](@ref), scaled up as
the total number of cells `n_cells` this frame falls short of `n_reference`
(the final, densest frame's cell count) -- `base_thresh` already tuned to
look good at `n_reference`. Uses a gentle `n^(-1/5)` power law (the same
scaling Silverman's rule of thumb uses for kernel-density bandwidth
selection, chosen for the same reason: density-based smoothing should widen
slowly with sparsity, not linearly or as `1/sqrt(n)`, which overshoots badly
at very low `n` -- e.g. only 4 cells at t=0). Confirmed numerically this
keeps t=180-360 (hundreds of cells) close to the already-tuned 40-45 while
softening t=0-90 (4-53 cells), where a flat 40 was comparatively too narrow
(that frame's own cells are much sparser, so the same absolute width
captures far fewer of them) and gave a jagged profile (max sample-to-sample
jump ~15.5) that adaptive scaling reduces to ~7.
"""
adaptive_ap_plane_thresh(n_cells, n_reference; base_thresh=40.0, power=0.2) =
    base_thresh * (n_reference / max(n_cells, 1))^power

"""
    smooth_radii(vecs::Vector{Vec3f}; window=9)

Smooth the MAGNITUDE of each vector in `vecs` with a centered moving average
over `window` samples (truncated at the ends, so the first/last points
average over fewer neighbors rather than wrapping or padding with zeros),
while leaving each vector's DIRECTION untouched -- i.e. this smooths the
radius profile specifically, not the left-right orientation. Even after
wider and adaptively-widened slices (see [`points_near_ap_plane`](@ref),
[`adaptive_ap_plane_thresh`](@ref)), the 100 independently-measured local
radii can still show sample-to-sample noise; this is an explicit smoothing
pass on the resulting profile itself, on top of (not instead of) that
spatial averaging.
"""
function smooth_radii(vecs::Vector{Vec3f}; window::Int=9)
    n = length(vecs)
    radii = norm.(vecs)
    # Guard against a degenerate (near-zero-radius) vector at an exact tip --
    # normalize would otherwise produce NaN, which propagates into the whole
    # smoothed sequence via the sum below.
    directions = [r > 1f-6 ? v / r : Vec3f(0, 0, 0) for (v, r) in zip(vecs, radii)]
    half = window ÷ 2
    smoothed = map(1:n) do j
        lo = max(1, j - half)
        hi = min(n, j + half)
        sum(@view radii[lo:hi]) / (hi - lo + 1)
    end
    return [s * d for (s, d) in zip(smoothed, directions)]
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

# Reference cell count for adaptive_ap_plane_thresh: the final (densest)
# frame's total cell count, since that's what the flat thresh=40 default was
# tuned against.
const N_CELLS_REFERENCE = length(get_pretwitch_points_at_time(pretwitch_df, series.frames[end].time))

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
    # between an EARLY vector and a LATE vector at each of the 100 sample
    # positions.
    #
    # EARLY: half the left-right extent measured at a LOCAL cross-sectional
    # slice through that sample's own AP position (not one fixed plane, and
    # NOT the widest point anywhere along the whole body) -- so the resulting
    # cylinder genuinely tapers near the head/tail, matching the real ovoid
    # shape of the point cloud, instead of applying one constant radius
    # everywhere.
    #
    # LATE: the actual local left-right half-distance/direction at each
    # terminal seam-cell pair, from right_vector_spline.
    #
    # Blending the full vector (not radius and direction separately) keeps
    # this identical in spirit to how interp_curve itself blends full
    # positions.
    adaptive_thresh = adaptive_ap_plane_thresh(length(all_pts), N_CELLS_REFERENCE)
    early_right_vecs = map(range(0, 1, length=100)) do x
        ap_position = lo + (hi - lo) * x
        slice_pts = points_near_ap_plane(all_pts, ref, ap_position; initial_thresh=adaptive_thresh)
        lr_projections = [dot(p - ref.origin, ref.lr) for p in slice_pts]
        lr_lo, lr_hi = extrema(lr_projections)
        Vec3f(((lr_hi - lr_lo) / 2) * ref.lr)
    end
    right_vecs = if isnothing(frame.right_vector_spline)
        early_right_vecs
    else
        [Vec3f((1 - w) * erv + w * Vec3f(frame.right_vector_spline(x)))
         for (erv, x) in zip(early_right_vecs, range(0, 1, length=100))]
    end
    surface = tube_mesh(interp_curve, smooth_radii(right_vecs))

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
