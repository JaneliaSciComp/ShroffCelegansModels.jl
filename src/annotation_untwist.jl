"""
    nearest_central_pt(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, expansion_factor = 1.0)

Find the nearest central point on the model to each annotation point.

# Arguments
- `model::AbstractCelegansModel`: The model to use for finding the nearest central point.
- `pts::AbstractVector{<: Point}`: The annotation points to find the nearest central point.
- `expansion_factor::Union{Vector{<: Real}, Real}=1.0`: The expansion factor to use for the distance check.
- `fourier_contour::Bool=false`: If `true`, estimate the outer contour of the cross section at each
  point's central angle by Fourier interpolation of the radial segments (see
  [`FourierContour`](@ref)). Two-step bound for speed: points beyond the sampled maximum radius
  are checked against the (lazily cached) interpolated maximum, and only points within it
  are compared with the contour at their own angle.

# Returns
- A tuple of three vectors: the spline parameter, the nearest central point, and the distance to the nearest central point.
"""
function nearest_central_pt(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, expansion_factor = 1.0; fourier_contour::Bool = false)
    # Sample points along the central spline
    cs = ShroffCelegansModels.central_spline(model)
    Npts = length(model)
    z = LinRange(0, 1, Npts)
    central_pts = cs.(z)

    # For each sampled point calculate the maximum radius of the cross section
    max_r_func = ShroffCelegansModels.max_radius_function(model)
    max_radius = max_r_func.(z)
    contour = fourier_contour ? FourierContour(model, z) : nothing

    # For each annotation point, find the nearest central point
    arg_pt_dist = map(pts) do pt
        # Distance between the annotation point and each sampled point along the central spline
        unbounded_dist = norm.(central_pts .- pt)
        a = 1
        if isnothing(contour)
            dist = copy(unbounded_dist)
            for ef in expansion_factor
                # Set all distances larger than the radius times the expansion factor to infinity
                dist[unbounded_dist .> max_radius .* ef] .= Inf
                # Find the index of the central point with the smallest distance
                a = argmin(dist)
                if dist[a] != Inf
                    # If we found a point, then we are done with this point
                    break
                end
                # Otherwise, reset the distances and try with a larger expansion factor
                dist = copy(unbounded_dist)
            end
            z[a], central_pts[a], dist[a]
        else
            best_i, best_d = 1, Inf
            for ef in expansion_factor
                best_i, best_d = 1, Inf
                for i in eachindex(z)
                    d = unbounded_dist[i]
                    # Cannot beat the best candidate so far
                    d < best_d || continue
                    # Step 1: beyond the maximum of the sampled radial segments?
                    # Then the true (interpolated) maximum is needed, computed lazily and cached.
                    if d > max_radius[i] * ef
                        # Step 2: still outside the true maximum: no further evaluation
                        d > true_max_radius(contour, i) * ef && continue
                    end
                    # Step 3: compare against the contour in the direction of the point
                    d <= contour_radius(contour, i, pt) * ef || continue
                    best_i, best_d = i, d
                end
                isfinite(best_d) && break
            end
            # Nothing within any bound: fall back to the unbounded nearest point like the original
            if !isfinite(best_d)
                best_i = argmin(unbounded_dist)
                best_d = Inf
            end
            z[best_i], central_pts[best_i], best_d
        end
    end
    second = x->x[2]
    # Split each tuple into three vectors
    return first.(arg_pt_dist), second.(arg_pt_dist), last.(arg_pt_dist)
end

"""
    get_central_point_parameters(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, central_spline_voxel_distances::AbstractVector{<: Real})

Get the parameters of the central spline for each annotation point given a
central spline voxel distance.
"""
function get_central_point_parameters(
    model::AbstractCelegansModel,
    pts::AbstractVector{<: Point},
    central_spline_voxel_distances::AbstractVector{<: Real}
)
    cs = ShroffCelegansModels.central_spline(model)
    Npts = length(model)
    z = LinRange(0, 1, Npts)
    central_pts = cs.(z)
    central_pts_dist = [0; cumsum(norm.(diff(central_pts)))]
    indices = map(central_spline_voxel_distances) do d
        findlast(<=(d + 0.1), central_pts_dist)
    end
    spline_param = z[indices]
    central_pts = central_pts[indices]
    pt_to_central_dist = norm.(central_pts .- pts)
    return spline_param, central_pts, pt_to_central_dist
end

"""
    FourierContour(model::AbstractCelegansModel, z)

Fourier-interpolated estimate of the outer contour of the worm at the spline
parameters `z`. At each `z` the lengths of the radial segments (transverse
splines) are treated as uniformly spaced samples of a periodic function of the
central angle. The angle of a point is measured from the first transverse
spline about the central spline tangent, the same convention used by
[`untwist_annotations`](@ref).

Use [`contour_radius`](@ref) for the radius in the direction of a point and
[`true_max_radius`](@ref) for the interpolated maximum, which is computed
lazily and cached per `z` index.
"""
struct FourierContour
    central_pts::Vector{Point3d}
    tangents::Vector{Point3d}
    rights::Vector{Point3d}
    coeffs::Matrix{ComplexF64}   # Fourier coefficients, N x length(z)
    freqs::Vector{Int}
    weights::Vector{Float64}
    orient::Float64
    max_cache::Vector{Float64}   # NaN until computed
end

_contour_angle(v, r, t) = atan(t ⋅ (r × v), r ⋅ v)

function FourierContour(model::AbstractCelegansModel, z)
    cs = central_spline(model)
    dcs = Derivative(1) * cs
    tss = transverse_splines(model)
    N = length(tss)
    central_pts = Point3d.(cs.(z))
    tangents = [Point3d(normalize(dcs(zi))) for zi in z]
    rights = [Point3d(normalize(tss[1](zi) - c)) for (zi, c) in zip(z, central_pts)]
    radii = [norm(tss[k](zi) - c) for k in 1:N, (zi, c) in zip(z, central_pts)]
    # Orientation of increasing spline index relative to the angle convention
    zmid = length(z) ÷ 2 + 1
    v2 = N > 1 ? tss[2](z[zmid]) - central_pts[zmid] : rights[zmid]
    o = sign(_contour_angle(v2, rights[zmid], tangents[zmid]))
    orient = o == 0 ? 1.0 : Float64(o)
    freqs = [m <= N ÷ 2 ? m : m - N for m in 0:N-1]
    # Split the Nyquist term for even N so the series stays real
    weights = [(iseven(N) && abs(m) == N ÷ 2) ? 0.5 : 1.0 for m in freqs]
    return FourierContour(central_pts, tangents, rights, fft(radii, 1) ./ N,
        freqs, weights, orient, fill(NaN, length(z)))
end

@inline function _contour_eval(c::FourierContour, i::Int, θ::Real)
    acc = 0.0
    @inbounds for j in eachindex(c.freqs)
        acc += c.weights[j] * real(c.coeffs[j, i] * cis(c.freqs[j] * θ))
    end
    return acc
end

"""
    contour_radius(c::FourierContour, i::Integer, pt) -> Float64

Interpolated contour radius at `z` index `i` in the direction of `pt`.
"""
function contour_radius(c::FourierContour, i::Integer, pt)
    θ = c.orient * _contour_angle(pt - c.central_pts[i], c.rights[i], c.tangents[i])
    return _contour_eval(c, i, θ)
end

"""
    true_max_radius(c::FourierContour, i::Integer) -> Float64

Maximum of the interpolated contour over all central angles at `z` index `i`,
found by dense sampling and local refinement. Cached after the first call.
"""
function true_max_radius(c::FourierContour, i::Integer)
    v = c.max_cache[i]
    isnan(v) || return v
    N = length(c.freqs)
    M = 16 * N
    step = 2π / M
    best_θ, best = 0.0, -Inf
    for k in 0:M-1
        r = _contour_eval(c, i, k * step)
        if r > best
            best, best_θ = r, k * step
        end
    end
    # Golden-section refinement within one grid step of the best sample
    lo, hi = best_θ - step, best_θ + step
    φ = (sqrt(5) - 1) / 2
    for _ in 1:30
        x1 = hi - φ * (hi - lo)
        x2 = lo + φ * (hi - lo)
        if _contour_eval(c, i, x1) > _contour_eval(c, i, x2)
            hi = x2
        else
            lo = x1
        end
    end
    best = max(best, _contour_eval(c, i, (lo + hi) / 2))
    c.max_cache[i] = best
    return best
end

function max_radius_function(model)
    function max_radius(z)
        tss = transverse_splines(model)
        central_pt = central_spline(model)(z)
        radii = map(tss) do ts
            norm(ts(z) - central_pt)
        end 
        return maximum(radii)
    end
end

function nearest_central_plane(model::AbstractCelegansModel, pts::AbstractVector{<: Point})
    cs = ShroffCelegansModels.central_spline(model)
    dcs = Derivative(1)*cs

    Npts = length(model)
    r = LinRange(0, 1, Npts)
    central_pts = cs.(r)


    ts1 = ShroffCelegansModels.transverse_spline(model, 1)
    right = ts1.(r)
    right_vec = right .- ncp
    right_norm = norm.(right_vec)
    # right_vec_unit = right_vec ./ right_norm
    # TODO: change to maximum radial distance
    max_radial_dist = right_norm

    normal_vecs = normalize.(dcs.(r))

    d = central_pts .⋅ normal_vecs

    arg_pt_dist = map(pts) do pt
        dist_to_plane = abs.(d .- (normal_vecs .⋅ pt))
        dist_to_central = norm.(central_pts .- pt)
        dist_along_plane = dist_to_central.^2 - dist_to_plane.^2
        a = argmin(dist_to_central)
        r[a], central_pts[a], dist_to_plane[a]
    end
    second = x->x[2]
    first.(arg_pt_dist), second.(arg_pt_dist), last.(arg_pt_dist)
end

"""
    untwist_annotations(model::AbstractCelegansModel, pts::AbstractVector{<: Point})

Untwist the annotations by computing their position relative to the central spline of the model.

# Arguments

- `model::AbstractCelegansModel`: The model to use for untwisting.
- `pts::AbstractVector{<: Point}`: The annotation points to untwist.
- `fourier_contour::Bool=false`: Use Fourier interpolation of the radial segments over the
  central angle to estimate the outer contour that bounds the anterior-posterior assignment,
  instead of the maximum radial segment length.

# Returns
- `Vector{Point3}`: The untwisted annotation points.

# Example
```julia
using ShroffCelegansModels
using ShroffCelegansModels.MIPAVIO

dataset = ShroffCelegansModels.NormalizedDataset("X:/shrofflab/OD1599_NU/120619_Pos2/Decon_reg/RegB")
mts = ShroffCelegansModels.ModelTimeSeries(dataset)
model = mts(1)
df = MIPAVIO.get_integrated_annotations(dataset, 1)
pts = MIPAVIO.mipav_df_to_points(df)
ShroffCelegansModels.untwist_annotations(model, pts)
```

# Notes
- The function uses the central spline of the model to untwist the annotations.
- The function uses the first transverse spline of the model to determine the right direction.
- The function uses the derivative of the central spline to determine the normal direction.

# Algorithm
1. Compute the nearest point on the central spline to each annotation point.
2. Compute the angle between the right direction and the annotation point.
3. Compute the distance between the central spline and the annotation point.
4. Compute the untwisted annotation point using the angle and distance.
5. Return the untwisted annotation points.
"""
function untwist_annotations(
    model::AbstractCelegansModel,
    pts::AbstractVector{<: Point},
    central_spline_voxel_distances::Union{AbstractVector{<: Real},Nothing} = nothing;
    fourier_contour::Bool = false
)::Vector{Point3d}
    thresholds = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.5, 2.0, 2.5]

    # For each point, find the nearest point along the central spline
    # t is the spline parameter
    # ncp is the nearest central point
    # pts_norm is the distance from the point to the central spline
    if isnothing(central_spline_voxel_distances)
        t, ncp, pts_norm = nearest_central_pt(model, pts, thresholds; fourier_contour)
    else
        t, ncp, pts_norm = get_central_point_parameters(model, pts, central_spline_voxel_distances)
    end

    # Objective: Compute the angle between the right direction and the point
    # The first transverse spline is used to determine which direction is the right side
    ts1 = ShroffCelegansModels.transverse_spline(model, 1)
    right = ts1.(t)

    # Obtain the central spline and its derivative
    # The derivative should point along the spline
    cs = ShroffCelegansModels.central_spline(model)
    dcs = Derivative(1)*cs

    # Obtain the vector from the nearest central point to the annotation
    pts_vec = pts .- ncp
    pts_vec_unit = pts_vec ./ pts_norm

    # Obtain the vector from the nearest central point to the right direction
    right_vec = right .- ncp
    right_norm = norm.(right_vec)
    right_vec_unit = right_vec ./ right_norm

    # Compute the cross product between the right direction and the annotation direction
    c = right_vec_unit .× pts_vec_unit
    # Compute the sign of the cross product and the derivative of the central spline
    # This indicates whether the cross product is pointing 
    s = sign.(dcs.(t) .⋅ c)

    # Compute the angle between the right direction and the annotation direction
    # by taking the arctangent of signed magnitude of the cross product
    # with the dot product of the right direction and the annotation direction
    angles = atan.(norm.(c) .*s, right_vec_unit .⋅ pts_vec_unit)

    # Obtain the straightened, untwisted, model
    smodel = StraightenedCelegansModel(model)
    scs = central_spline(smodel)
    # Compute the z-coordinate of the nearest central point on the straightened model
    z = last.(scs.(t))

    # The returned untwisted annotation points are computed using the angle and distance
    # between the central spline and the annotation point
    # The z-coordinate is the same as that of the nearest central point on the straightened model
    return map(zip(angles,pts_norm,z)) do (angle, dist, z)
        Point3d(cos(angle)*dist, sin(angle)*dist, z)
    end
end

function untwist_annotation(
    model::AbstractCelegansModel,
    pt::Point,
    central_spline_voxel_distance::Union{Real,Nothing} = nothing # "z"
)
    pts = [pt]
    central_spline_voxel_distances =
        isnothing(central_spline_voxel_distance) ?
            nothing :
            [central_spline_voxel_distance]
    return only(untwist_annotations(model, pts, central_spline_voxel_distances))
end

"""
    untwist_annotations(dataset::NormalizedDataset, timepoint::Int=1)

Untwist the annotations at a specific timepoint in the dataset.

# Arguments
- `dataset::NormalizedDataset`: The dataset to use for untwisting.
- `timepoint::Int=1`: The timepoint to use for untwisting, where 1 indicates of the first time point of the dataset.

# Returns
- `Dict{String, Vector{Point3}}`: The untwisted annotations using the annotation name as a the key

# Example
```julia
using ShroffCelegansModels

dataset = ShroffCelegansModels.NormalizedDataset("X:/shrofflab/OD1599_NU/120619_Pos2/Decon_reg/RegB")
# Obtain the untwisted annotations for the C3 cell at timepoint 71
c3_pt = ShroffCelegansModels.untwist_annotations(dataset, 71)["C3"]
```
"""
function untwist_annotations(dataset::NormalizedDataset, timepoint::Int=1)::Union{Missing, Dict{String, Point3d}}
    mts = ModelTimeSeries(dataset)
    try
        df = MIPAVIO.get_integrated_annotations(dataset, timepoint)
        pts = MIPAVIO.mipav_df_to_points(df)::Vector{Point3d}
        return Dict{String,Point3d}(Vector{String}(df.name)::Vector{String} .=> untwist_annotations(mts(timepoint), pts))
    catch err
        return missing
    end
end

function twisted_annotations(dataset::NormalizedDataset, timepoint::Int=1)
    try
        df = MIPAVIO.get_integrated_annotations(dataset, timepoint)
        pts = MIPAVIO.mipav_df_to_points(df)
        return Dict(df.name .=> pts)
    catch err
        return missing
    end
end

function distance_to_twisted_annotation(model::AbstractCelegansModel, pt::Point)
    # Sample points along the central spline
    cs = ShroffCelegansModels.central_spline(model)
    dcs = Derivative(1)*cs
    Npts = length(model)
    z = LinRange(0, 1, Npts)
    central_pts = cs.(z)
    annotation_vecs = -(central_pts .- pt)
    annotation_norms = norm.(annotation_vecs)
    unit_annotation_vecs = annotation_vecs ./ annotation_norms

    ts1 = ShroffCelegansModels.transverse_spline(model, 1)
    right = ts1.(z)
    right_vecs = right .- central_pts
    right_norms = norm.(right_vecs)
    unit_right_vecs = right_vecs ./ right_norms

    # Compute the cross product between the right direction and the annotation direction
    c = unit_right_vecs .× unit_annotation_vecs

    central_spline_vecs = dcs.(z)
    central_spline_norms = norm.(central_spline_vecs)
    unit_central_spline_vecs = central_spline_vecs ./ central_spline_norms

    # Compute the sign of the cross product and the derivative of the central spline
    # This indicates whether the cross product is pointing 
    s = sign.(unit_central_spline_vecs .⋅ c)

    # by taking the arctangent of signed magnitude of the cross product
    # with the dot product of the right direction and the annotation direction
    angles = atan.(norm.(c) .*s, unit_right_vecs .⋅ unit_annotation_vecs)

    return (; annotation_norms, angles, unit_annotation_vecs, unit_right_vecs, unit_central_spline_vecs, c, central_pts, right_vecs, right_norms)
end