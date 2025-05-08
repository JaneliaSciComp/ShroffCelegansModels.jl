"""
    nearest_central_pt(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, expansion_factor = 1.0)

Find the nearest central point on the model to each annotation point.

# Arguments
- `model::AbstractCelegansModel`: The model to use for finding the nearest central point.
- `pts::AbstractVector{<: Point}`: The annotation points to find the nearest central point.
- `expansion_factor::Union{Vector{<: Real}, Real}=1.0`: The expansion factor to use for the distance check.

# Returns
- A tuple of three vectors: the spline parameter, the nearest central point, and the distance to the nearest central point.
"""
function nearest_central_pt(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, expansion_factor = 1.0)
    # Sample points along the central spline
    cs = ShroffCelegansModels.central_spline(model)
    Npts = length(model)
    z = LinRange(0, 1, Npts)
    central_pts = cs.(z)

    # For each sampled point calculate the maximum radius of the cross section
    max_r_func = ShroffCelegansModels.max_radius_function(model)
    radius = max_r_func.(z)

    # For each annotation point, find the nearest central point
    arg_pt_dist = map(pts) do pt
        # Distance between the annotation point and each sampled point along the central spline
        unbounded_dist = norm.(central_pts .- pt)
        # bound check
        a = 1
        dist = copy(unbounded_dist)
        for ef in expansion_factor
            # Set all distances larger than the radius times the expansion factor to infinity
            dist[unbounded_dist .> radius .* ef] .= Inf
            # Find the index of the central point with the smallest distance
            a = argmin(dist)
            if dist[a] != Inf
                # If we found a point, then we are done with this point
                break
            end
            # Otherwise, reset the distances and try with a larger expansion factor
            dist = copy(unbounded_dist)
        end
        # Return a tuple of the spline parameter, the nearest central point, and the distance to the nearest central point
        z[a], central_pts[a], dist[a]
    end
    second = x->x[2]
    # Split each tuple into three vectors
    return first.(arg_pt_dist), second.(arg_pt_dist), last.(arg_pt_dist)
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
function untwist_annotations(model::AbstractCelegansModel, pts::AbstractVector{<: Point})
    thresholds = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.5, 2.0, 2.5]

    # For each point, find the nearest point along the central spline
    # t is the spline parameter
    # ncp is the nearest central point
    # pts_norm is the distance from the point to the central spline
    t, ncp, pts_norm = nearest_central_pt(model, pts, thresholds)

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
    return map(zip(angles,pts_norm,z)) do (angle, dist,z)
        Point3(cos(angle)*dist, sin(angle)*dist, z)
    end
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
function untwist_annotations(dataset::NormalizedDataset, timepoint::Int=1)
    mts = ModelTimeSeries(dataset)
    try
        df = MIPAVIO.get_integrated_annotations(dataset, timepoint)
        pts = MIPAVIO.mipav_df_to_points(df)
        return Dict(df.name .=> untwist_annotations(mts(timepoint), pts))
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