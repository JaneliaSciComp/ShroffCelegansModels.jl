function cross_sections_at_knots(model)
    pts = interpolation_points(model.central_spline)
    sections = map(transverse_splines(model)) do spline
        cs_pts = spline.(pts)
    end
    sections = stack(sections)
    sections = [sections @view(sections[:,1]) fill(Point3(NaN), size(sections,1), 1)]
    sections = PermutedDimsArray(sections, (2,1))
    return vec(sections)
end
