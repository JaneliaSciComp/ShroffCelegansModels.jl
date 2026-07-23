# Establish anterior-posterior (AP), dorsal-ventral (DV), and left-right (LR)
# axes for the pretwitch embryo (see pretwitch_orientation.txt for the design
# notes this implements).
#
# Approach 1 (WormGuides compass interpolation): the AP axis is fixed, derived
# once from the four-cell stage (ABa/ABp/P2/EMS). The DV/LR basis is assumed
# to rotate rigidly about that fixed AP axis over the course of development,
# with the rotation angle given by WormGuides' Key_Frames_Rotate /
# Key_Values_Rotate compass table. The compass's 90 degree value at frame 1
# is defined to reproduce the directly observed four-cell-stage DV direction,
# so the interpolated angle is a rotation *relative to that calibration*, not
# an absolute lab-frame heading.
#
# This is explicitly a first, unvalidated attempt (see pretwitch_orientation.txt
# item 1) -- the left/right sign of the LR axis is NOT anatomically resolved
# here (that requires either lineage-descendant markers, approach 2, or
# cross-checking against posttwitch orientation as in lattice_orientation.jl,
# approach 3).

const _WORMGUIDES_COMPASS_KEY_FRAMES = [1, 16, 321, 359]
const _WORMGUIDES_COMPASS_KEY_VALUES_DEG = [90.0, 30.0, 30.0, 90.0]

"""
    wormguides_compass_angle(frame; key_frames=_WORMGUIDES_COMPASS_KEY_FRAMES,
                                     key_values_deg=_WORMGUIDES_COMPASS_KEY_VALUES_DEG)

Piecewise-linear interpolation of the WormGuides compass rotation angle (in
radians) at a given `frame`, per the `Key_Frames_Rotate` / `Key_Values_Rotate`
table in pretwitch_orientation.txt. Clamped to the first/last key value outside
the key frame range.
"""
function wormguides_compass_angle(
    frame::Real;
    key_frames::AbstractVector{<:Real}=_WORMGUIDES_COMPASS_KEY_FRAMES,
    key_values_deg::AbstractVector{<:Real}=_WORMGUIDES_COMPASS_KEY_VALUES_DEG,
)
    if frame <= key_frames[1]
        return deg2rad(key_values_deg[1])
    elseif frame >= key_frames[end]
        return deg2rad(key_values_deg[end])
    end
    i = searchsortedlast(key_frames, frame)
    f0, f1 = key_frames[i], key_frames[i+1]
    v0, v1 = key_values_deg[i], key_values_deg[i+1]
    L = (frame - f0) / (f1 - f0)
    return deg2rad(v0 + L * (v1 - v0))
end

"""
    pretwitch_time_to_wormguides_frame(time)

Pretwitch `time` runs 0-360 (361 timepoints, matching `get_pretwitch_df`'s
`time` column); WormGuides frames run 1-360. Frame = time + 1.
"""
pretwitch_time_to_wormguides_frame(time::Real) = time + 1

"""
    pretwitch_four_cell_points(pretwitch_df; time=0)

Return the (ABa, ABp, P2, EMS) points at the four-cell stage (`time`, default
0) as a NamedTuple of `Point3f`s.
"""
function pretwitch_four_cell_points(pretwitch_df; time::Int=0)
    pts = get_pretwitch_points_at_time(pretwitch_df, time)
    return (ABa=pts["ABa"], ABp=pts["ABp"], P2=pts["P2"], EMS=pts["EMS"])
end

"""
    pretwitch_cells_by_name(pretwitch_df)

Group `pretwitch_df` by cell/lineage name, returning
`Dict{String, Vector{Tuple{Int,Point3f}}}` of that cell's own
`(time, position)` pairs (its existence window, before it either divides into
two differently-named daughters or the pretwitch dataset ends).
"""
function pretwitch_cells_by_name(pretwitch_df)
    grouped = Dict{String, Vector{Tuple{Int,Point3f}}}()
    for row in eachrow(pretwitch_df)
        push!(get!(() -> Tuple{Int,Point3f}[], grouped, row.cell),
              (row.time, Point3f(row.x, row.y, row.z)))
    end
    return grouped
end

"""
    pretwitch_lineage_track(target_name, grouped)

Reconstruct a single continuous position trajectory for `target_name` by
splicing together the existence windows of its full ancestor chain -- e.g.
for `"Cpaaaa"`, the chain is `"C"`, `"Cp"`, `"Cpa"`, `"Cpaa"`, `"Cpaaa"`,
`"Cpaaaa"` (each a division away from the next). `grouped` is the output of
[`pretwitch_cells_by_name`](@ref). Ancestor prefixes not present as their own
distinct name in `grouped` (division doesn't always add exactly one
character in a way that leaves every prefix separately observed) are simply
skipped. Returns `(times, points)` sorted by time.

This is what "working back from" a cell like Cpaaaa (identified in the
posttwitch lattice orientation QC, see `lattice_orientation.jl`) to the
pretwitch stage means concretely: its full ancestry, not just the single
frame where `target_name` itself first appears.
"""
function pretwitch_lineage_track(target_name::AbstractString, grouped::AbstractDict{String})
    track = Tuple{Int,Point3f}[]
    for k in 1:length(target_name)
        prefix = target_name[1:k]
        haskey(grouped, prefix) && append!(track, grouped[prefix])
    end
    sort!(track, by=first)
    return first.(track), last.(track)
end

"""
    pretwitch_founder_group(name::AbstractString)

Classify a pretwitch cell/lineage `name` by which of the four founder cells
(ABa, ABp, EMS, P2) it descends from, returning `:ABa`, `:ABp`, `:EMS`, or
`:P2`.

Unlike the AB lineage (which keeps growing an "AB" + a/p/l/r/d/v prefix
indefinitely, e.g. `ABplaaappa` -- see [`get_full_lineage_df`](@ref)), the
other two four-cell-stage founders are *renamed* at each division rather than
suffixed: EMS -> MS + E, and P2 -> C + P3 -> D + P4 -> Z2 + Z3. So membership
is determined by a fixed table of lineage-name prefixes rather than a single
common root string. Verified to classify all 1332 distinct cell names in
`final_20251205_pre-twitch_coords.csv` with no unknowns.
"""
function pretwitch_founder_group(name::AbstractString)
    startswith(name, "ABa") && return :ABa
    startswith(name, "ABp") && return :ABp
    (name == "EMS" || startswith(name, "E") || startswith(name, "MS")) && return :EMS
    (name == "P2" || startswith(name, "P3") || startswith(name, "P4") ||
     startswith(name, "Z2") || startswith(name, "Z3") ||
     startswith(name, "C") || startswith(name, "D")) && return :P2
    error("pretwitch_founder_group: unrecognized cell name $(repr(name))")
end

"""
    PRETWITCH_FOUNDER_COLORS

Default display color for each of the four founder-cell lineage groups (see
[`pretwitch_founder_group`](@ref)), chosen to be visually distinct from the
red/blue/green used elsewhere for the AP/DV/LR axes themselves.
"""
const PRETWITCH_FOUNDER_COLORS = Dict(
    :ABa => :orange,
    :ABp => :purple,
    :EMS => :cyan,
    :P2 => :magenta,
)

"""
    PretwitchReferenceAxes

Orthonormal (`ap`, `dv`, `lr`) unit vectors and the four-cell-stage `origin`
(centroid of ABa/ABp/P2/EMS) that anchor the pretwitch orientation scheme.
`ap` points anterior -> posterior, `dv` points ventral -> dorsal (after
orthogonalizing against `ap`), and `lr` completes a right-handed frame via
`cross(ap, dv)` -- its correspondence to anatomical left/right is NOT
resolved by this calibration alone.
"""
struct PretwitchReferenceAxes
    origin::Point3f
    ap::Vec3f
    dv::Vec3f
    lr::Vec3f
end

"""
    _pretwitch_axes_from_four_points(ABa, ABp, P2, EMS)

Shared orthonormal-basis construction used by both the four-cell-stage
calibration ([`pretwitch_reference_axes`](@ref), Approach 1) and the
per-frame lineage-midpoint scheme ([`pretwitch_axes_from_group_midpoints`](@ref),
Approach 2): `ap` points anterior (ABa) -> posterior (P2), `dv` points ventral
(EMS) -> dorsal (ABp) after orthogonalizing against `ap`, and `lr` completes a
right-handed frame via `cross(ap, dv)`.
"""
function _pretwitch_axes_from_four_points(ABa, ABp, P2, EMS)
    origin = Point3f((ABa .+ ABp .+ P2 .+ EMS) ./ 4)

    ap = normalize(P2 - ABa)
    dv_raw = ABp - EMS
    dv = normalize(dv_raw - dot(dv_raw, ap) * ap)
    lr = normalize(cross(ap, dv))

    return PretwitchReferenceAxes(origin, ap, dv, lr)
end

function pretwitch_reference_axes(pretwitch_df; time::Int=0)
    cells = pretwitch_four_cell_points(pretwitch_df; time)
    return _pretwitch_axes_from_four_points(cells.ABa, cells.ABp, cells.P2, cells.EMS)
end

"""
    pretwitch_axes_from_group_midpoints(group_mid)

Approach 2 (see pretwitch_orientation.txt item 2): instead of interpolating a
fixed-AP-axis compass angle (Approach 1), derive the AP/DV/LR axes fresh at
*every* frame from the current centroid ("midpoint") of each founder's
descendants (see [`pretwitch_founder_group`](@ref)). `group_mid` maps
`:ABa`/`:ABp`/`:EMS`/`:P2` to that group's current-frame centroid `Point3f`
(the same centroids already used for the star markers in the movie).

Unlike Approach 1, the AP axis here is NOT fixed -- it tracks however the
ABa/P2 descendant midpoints actually move, so this scheme captures any
translation/curvature of the true anterior-posterior axis that the fixed-AP
compass interpolation cannot.
"""
function pretwitch_axes_from_group_midpoints(group_mid::AbstractDict{Symbol})
    return _pretwitch_axes_from_four_points(
        group_mid[:ABa], group_mid[:ABp], group_mid[:P2], group_mid[:EMS])
end

"""
    PretwitchAxesAtTime

`ap`/`dv`/`lr` unit vectors for a single pretwitch timepoint, with `dv`/`lr`
rotated about the fixed `reference.ap` axis by the interpolated WormGuides
compass angle for that timepoint (see [`wormguides_compass_angle`](@ref)).
"""
struct PretwitchAxesAtTime
    time::Int
    ap::Vec3f
    dv::Vec3f
    lr::Vec3f
    compass_angle_rad::Float64
end

"""
    pretwitch_axes_at_time(reference::PretwitchReferenceAxes, time)

Rotate `reference`'s DV/LR basis about the fixed AP axis by the WormGuides
compass angle interpolated for `time`. At `time` equal to the reference's own
four-cell-stage timepoint (compass angle 90 degrees), this reproduces
`reference.dv`/`reference.lr` exactly by construction.
"""
function pretwitch_axes_at_time(reference::PretwitchReferenceAxes, time::Int)
    frame = pretwitch_time_to_wormguides_frame(time)
    θ = wormguides_compass_angle(frame)
    # vector(θ) = cos(θ) * lr + sin(θ) * dv, so θ=90deg reproduces `dv` and
    # calibrates the compass against the directly observed four-cell axes.
    dv = cos(θ) * reference.lr + sin(θ) * reference.dv
    lr = normalize(cross(reference.ap, dv))
    return PretwitchAxesAtTime(time, reference.ap, dv, lr, θ)
end
