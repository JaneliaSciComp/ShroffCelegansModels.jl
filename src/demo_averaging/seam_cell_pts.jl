using BSplineKit.SplineInterpolations: interpolation_points

function seam_cell_pts(model, n_upsample)
    pts = interpolation_points(model)[1:2^n_upsample:end]
    _transverse_splines = transverse_splines(model)
    [_transverse_splines[1].(pts); _transverse_splines[17].(pts)]
end