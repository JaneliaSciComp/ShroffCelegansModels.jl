# Pretwitch central-spline construction (see pretwitch_orientation.txt and
# this session's design discussion for the full rationale).
#
# Mirrors src/build_model.jl's posttwitch central spline (fit through the
# midpoints of 10 left/right seam-cell pairs), but pretwitch data has no
# seam cells yet -- only their embryonic-lineage ancestors. Ancestor
# positions are reconstructed via pretwitch_lineage_track (splicing together
# the ancestor chain, as seam_cell_to_lineage_map.jl's
# get_full_lineage_seam_cell_points already does for a single seam cell).
#
# The spline is built ONLY from TERMINAL seam cells -- pairs whose ancestor
# tracking has reached the exact final named seam-cell identity (both L and
# R), not merely some shorter, still-dividing ancestor prefix. Connecting
# undifferentiated shared ancestors (e.g. H0 and H1 both still being nothing
# more than "whatever ABa/ABp average to" before either has committed to
# becoming a seam cell) doesn't represent real seam-cell geometry -- it's a
# stand-in for cells that don't exist as such yet. Once restricted to
# terminal cells, no two stations can ever coincide (each is a real, distinct
# differentiated cell), so no merge/dedup step is needed at all.
#
# The practical cost: terminal seam cells don't exist until quite late.
# Confirmed empirically: ZERO of the 10 pairs are terminal for frames 0-169
# (the first 47% of the pretwitch window); the count then climbs
# non-monotonically and only reaches its max (9 of 10 -- V5's pair never
# reaches its exact terminal name within this dataset, stopping one division
# short) from frame ~204 onward. So for roughly the first half of the movie
# there may be too few terminal knots (or none) to fit any spline at all --
# `PretwitchCentralSplineFrame.central_spline` is `nothing` in that case.

"""
    PretwitchPairTrack

Full per-frame (0:360) reconstruction for one canonical AP-order left/right
seam-cell pair (e.g. "H0"). `left_name`/`right_name` are the currently
active spliced ancestor name backing `left[t+1]`/`right[t+1]` (see
[`pretwitch_lineage_track`](@ref)); `left_name[i] == right_name[i]` means
this pair has not yet diverged into two distinct cells at that frame.
`target_left`/`target_right` are the exact final seam-cell lineage names
(from `seam_cell_to_lineage_map`) that `left_name[i]`/`right_name[i]` must
equal for this pair to count as TERMINAL (a real, fully differentiated seam
cell) at frame `i`.
"""
struct PretwitchPairTrack
    name::String
    left_name::Vector{String}
    right_name::Vector{String}
    left::Vector{Point3f}
    right::Vector{Point3f}
    target_left::String
    target_right::String
end

"""
    is_terminal(track::PretwitchPairTrack, i::Int)

Whether pair `track` is a real, fully differentiated (not merely a shared,
still-dividing ancestor) seam cell at frame index `i` -- i.e. both
`left_name[i]` and `right_name[i]` have reached their exact final lineage
name.
"""
is_terminal(track::PretwitchPairTrack, i::Int) =
    track.left_name[i] == track.target_left && track.right_name[i] == track.target_right

"""
    pretwitch_pair_tracks(pretwitch_df=get_pretwitch_df())

Build a [`PretwitchPairTrack`](@ref) for each of the 10 canonical L/R
seam-cell pairs (`left_seam_cells`/`right_seam_cells`, anterior->posterior
order), reconstructing each side's full ancestor-spliced position track via
[`pretwitch_lineage_track`](@ref).
"""
function pretwitch_pair_tracks(pretwitch_df=get_pretwitch_df())
    grouped = pretwitch_cells_by_name(pretwitch_df)
    map(zip(left_seam_cells, right_seam_cells)) do (l, r)
        target_left, target_right = seam_cell_to_lineage_map[l], seam_cell_to_lineage_map[r]
        lt, lp, ln = pretwitch_lineage_track(target_left, grouped)
        rt, rp, rn = pretwitch_lineage_track(target_right, grouped)
        @assert lt == rt "pretwitch_pair_tracks: L/R frame coverage mismatch for pair $l/$r"
        PretwitchPairTrack(replace(l, 'L' => ""), ln, rn, lp, rp, target_left, target_right)
    end
end

"""
    PairConfidence

Diagnostic-only confidence that a pair's midpoint approximates a true
bilateral average, for one pair at one frame. When `diverged` (left/right
resolve to different named cells), the midpoint is a real L/R average and is
always trusted (`score=1`). When not yet diverged, the midpoint is really
just that one shared, not-yet-split cell's own position -- `score` then
comes from how strongly the regression (`pretwitch_axis_regression_weights.csv`,
see scripts/regress_pretwitch_axis_weights.jl) found that specific cell's
position to correlate with the true LR axis: a cell already sitting well off
the midline (large `|lr_weight|`) is scored as less trustworthy as a midline
stand-in, even though structurally it looks "undiverged."
"""
struct PairConfidence
    diverged::Bool
    lr_weight_left::Float64
    lr_weight_right::Float64
    score::Float64
end

"""
    load_lr_weights(csv_path=<repo root>/pretwitch_axis_regression_weights.csv)

Load the `name => lr_weight` column of the regression CSV. Returns an empty
`Dict` (all names default to weight 0, i.e. `score=1`, fully trusted) if the
CSV doesn't exist -- it's gitignored/regenerated, not committed, so this
degrades gracefully on a fresh checkout instead of erroring.
"""
function load_lr_weights(csv_path=joinpath(@__DIR__, "..", "..", "pretwitch_axis_regression_weights.csv"))
    lr_w = Dict{String,Float64}()
    isfile(csv_path) || return lr_w
    open(csv_path) do io
        readline(io)  # header
        for line in eachline(io)
            parts = split(line, ',')
            lr_w[parts[1]] = parse(Float64, parts[4])
        end
    end
    return lr_w
end

function pair_confidence(track::PretwitchPairTrack, i::Int, lr_w::AbstractDict{String,Float64}, max_abs_weight::Float64)
    ln, rn = track.left_name[i], track.right_name[i]
    diverged = ln != rn
    wl, wr = get(lr_w, ln, 0.0), get(lr_w, rn, 0.0)
    score = diverged ? 1.0 : clamp(1.0 - abs(wl) / max_abs_weight, 0.0, 1.0)
    return PairConfidence(diverged, wl, wr, score)
end

"""
    PretwitchCentralSplineFrame

One pretwitch timepoint's central-spline state: the fitted natural spline
(parameterized by normalized arc length, `afTime ∈ [0,1]`, exactly as
`build_celegans_model`'s posttwitch `center_spline`) through whichever
stations are currently TERMINAL (see [`is_terminal`](@ref)), in fixed
NOMINAL fate order (H0,H1,H2,V1..V6,T -- never re-sorted by current spatial
position, so a knot's anatomical identity is consistent across every frame),
`terminal_stations` (1-based nominal indices of the knots actually used), and
per-pair diagnostic [`PairConfidence`](@ref) (still reported for all 10
pairs, independent of terminal status).

`central_spline` is `nothing` when fewer than 2 stations are terminal this
frame (no curve can be fit) -- true for frames 0-169 (~47% of the pretwitch
window), where ZERO pairs have differentiated into a real seam cell yet.

`is_ap_monotonic` is a diagnostic ONLY: whether this frame's nominal knot
order also happens to be spatially monotonic along the fixed AP reference
axis. It does not affect the fit. `false` means the seam-cell ancestors
haven't yet spatially sorted into their final AP body order at this
frame -- expected during gastrulation-stage rearrangement, not an error.
`missing` when fewer than 2 knots (monotonicity is undefined/vacuous).
"""
struct PretwitchCentralSplineFrame
    time::Int
    terminal_stations::Vector{Int}
    knot_positions::Vector{Point3f}
    afTime::Vector{Float64}
    central_spline::Union{Nothing,BSplineKit.SplineWrapper}
    pair_confidence::Vector{PairConfidence}
    is_ap_monotonic::Union{Missing,Bool}
end

"""
    _fit_central_spline(afTime, knot_positions)

Fit a natural spline through `knot_positions`, exactly as
[`ParametricSplines.interpolate_natural_cubic_spline`](@ref) does for the
posttwitch model (which always has exactly 10 knots) -- EXCEPT that a
terminal-only pretwitch frame can have as few as 2 knots, and `BSplineKit`'s
natural cubic (`BSplineOrder(4)`) interpolation is only evaluable with >=4
points (verified directly: it *constructs* with 3 but throws `ArgumentError:
wrong number of coefficients` when called). `Natural` boundary conditions are
also only supported for even spline orders. So: use cubic (order 4) when
there are enough knots, otherwise fall back to a natural linear spline
(order 2, i.e. a polyline/straight segment) -- a reasonable representation
for 2-3 widely-spaced knots anyway. Caller must ensure `length(knot_positions) >= 2`.
"""
function _fit_central_spline(afTime, knot_positions)
    order = length(knot_positions) >= 4 ? 4 : 2
    return BSplineKit.interpolate(afTime, knot_positions, BSplineKit.BSplineOrder(order), BSplineKit.Natural())
end

"""
    central_spline_frame(tracks, i, lr_w, max_abs_weight, ap_ref)

Build one [`PretwitchCentralSplineFrame`](@ref) at pair-track index `i`
(frame time `i-1`): select only currently-TERMINAL pairs (real, fully
differentiated seam cells -- see [`is_terminal`](@ref)), fit a natural spline
through their midpoints (see [`_fit_central_spline`](@ref)) when there are
enough of them, and compute per-pair diagnostic confidence for all 10 pairs
regardless of terminal status.
"""
function central_spline_frame(tracks::Vector{PretwitchPairTrack}, i::Int, lr_w::AbstractDict{String,Float64}, max_abs_weight::Float64, ap_ref::PretwitchReferenceAxes)
    n = length(tracks)
    terminal_stations = [k for k in 1:n if is_terminal(tracks[k], i)]
    knot_positions = [Point3f((tracks[k].left[i] .+ tracks[k].right[i]) ./ 2) for k in terminal_stations]

    if length(knot_positions) >= 2
        # Diagnostic only (does NOT affect the fit): does this frame's
        # nominal order also happen to be spatially monotonic along the
        # fixed AP reference axis? Early-to-mid development frequently
        # answers "no" -- the seam-cell ancestors are still undergoing
        # gastrulation-like rearrangement and are not yet sorted into their
        # final AP body order. That's a real developmental signal, not a bug
        # to fix by reordering.
        ap_projections = [dot(p - ap_ref.origin, ap_ref.ap) for p in knot_positions]
        is_ap_monotonic = issorted(ap_projections)

        deltaLengths = norm.(diff(knot_positions))
        afTime = Float64[0; cumsum(deltaLengths)]
        afTime ./= afTime[end]
        central_spline = _fit_central_spline(afTime, knot_positions)
    else
        is_ap_monotonic = missing
        afTime = Float64[]
        central_spline = nothing
    end

    confidences = [pair_confidence(tracks[k], i, lr_w, max_abs_weight) for k in 1:n]
    return PretwitchCentralSplineFrame(i - 1, terminal_stations, knot_positions, afTime, central_spline, confidences, is_ap_monotonic)
end

"""
    PretwitchCentralSplineSeries

`frames[i].time == i-1` for `i in 1:361`; `station_names` gives the 10
canonical AP-order pair names (`"H0"`, ..., `"T"`) that `pair_confidence`
(and `terminal_stations`, by 1-based index) refer to.
"""
struct PretwitchCentralSplineSeries
    frames::Vector{PretwitchCentralSplineFrame}
    station_names::Vector{String}
end

"""
    build_pretwitch_central_spline_series(pretwitch_df=get_pretwitch_df())

Build the full 361-frame [`PretwitchCentralSplineSeries`](@ref): per frame,
fit a central spline through whichever seam-cell pairs are currently
TERMINAL (real, fully differentiated cells -- see [`is_terminal`](@ref)).
No pair is terminal before frame 170; the count then climbs
non-monotonically, reaching 9 of 10 (V5 never reaches its exact terminal
name in this dataset) from around frame 204 onward. `central_spline` is
`nothing` for frames with fewer than 2 terminal pairs.
"""
function build_pretwitch_central_spline_series(pretwitch_df=get_pretwitch_df())
    tracks = pretwitch_pair_tracks(pretwitch_df)
    lr_w = load_lr_weights()
    max_abs_weight = isempty(lr_w) ? 1.0 : maximum(abs, values(lr_w))
    ap_ref = pretwitch_reference_axes(pretwitch_df)
    n_times = length(tracks[1].left)
    frames = [central_spline_frame(tracks, i, lr_w, max_abs_weight, ap_ref) for i in 1:n_times]
    return PretwitchCentralSplineSeries(frames, [tr.name for tr in tracks])
end
