function nearest_central_pt(model::AbstractCelegansModel, pts::AbstractVector{<: Point}, expansion_factor = 1.0)
    cs = ShroffCelegansModels.central_spline(model)
    Npts = length(model)
    z = LinRange(0, 1, Npts)
    central_pts = cs.(z)

    # bound the distance check

    #=
    ts1 = ShroffCelegansModels.transverse_spline(model, 1)
    ts1_pts = ts1.(r)
    radius1 = norm.(ts1_pts .- central_pts) .* expansion_factor
    =#
    max_r_func = ShroffCelegansModels.max_radius_function(model)
    radius = max_r_func.(z)

    arg_pt_dist = map(pts) do pt
        unbounded_dist = norm.(central_pts .- pt)
        # bound check
        a = 1
        dist = copy(unbounded_dist)
        for ef in expansion_factor
            dist[unbounded_dist .> radius .* ef] .= Inf
            a = argmin(dist)
            b = argmin(unbounded_dist)
            if dist[a] != Inf
                break
            end
            dist = copy(unbounded_dist)
        end
        z[a], central_pts[a], dist[a]
    end
    second = x->x[2]
    third = x->x[3]
    first.(arg_pt_dist), second.(arg_pt_dist), last.(arg_pt_dist)
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

function untwist_annotations(model::AbstractCelegansModel, pts::AbstractVector{<: Point})
    thresholds = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.5, 1.5, 2.0, 2.5]
    #thresholds = [1.0, 1.5, 1.5, 2.0, 2.5]
    t, ncp, pts_norm = nearest_central_pt(model, pts, thresholds)
    ts1 = ShroffCelegansModels.transverse_spline(model, 1)
    right = ts1.(t)

    cs = ShroffCelegansModels.central_spline(model)
    dcs = Derivative(1)*cs

    pts_vec = pts .- ncp
    pts_vec_unit = pts_vec ./ pts_norm

    right_vec = right .- ncp
    right_norm = norm.(right_vec)
    right_vec_unit = right_vec ./ right_norm

    c = right_vec_unit .× pts_vec_unit
    s = sign.(dcs.(t) .⋅ c)

    angles = atan.(norm.(c) .*s, right_vec_unit .⋅ pts_vec_unit)

    # return angles
    # return pts_norm

    # untwisted model
    smodel = StraightenedCelegansModel(model)
    scs = central_spline(smodel)
    z = last.(scs.(t))

    map(zip(angles,pts_norm,z)) do (angle, dist,z)
        Point3(cos(angle)*dist, sin(angle)*dist, z)
    end
end

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
