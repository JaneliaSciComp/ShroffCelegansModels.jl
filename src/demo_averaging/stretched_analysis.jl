function stretched_analysis(smodel)
    model = parent(smodel)
    r = LinRange(0,1,length(smodel))

    cs_pts = ShroffCelegansModels.central_spline(model).(r)
    scs_pts = ShroffCelegansModels.central_spline(smodel).(r)

    t1_pts = ShroffCelegansModels.transverse_spline(model, 1).(r)

    st1_pts = t1_pts .- cs_pts .+ scs_pts
    to_cylindrical = CoordinateTransformations.CylindricalFromCartesian()
    θ = to_cylindrical.(st1_pts) .|> x->x.θ
end