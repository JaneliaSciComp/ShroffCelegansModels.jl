function radial_cross_section(model::AbstractCelegansModel, p)
    _central_spline = central_spline(model)
    center = _central_spline(p)
    map(transverse_splines(model)) do transverse_spline
        norm(transverse_spline(p) - center)
    end
end

function cross_section_area(model::AbstractCelegansModel, p)
    radii = radial_cross_section(model, p)
    radii_hat = fft(radii)
    N = length(radii)
    _sum_of_square_norm = sum(radii_hat) do radius_hat
        abs2(radius_hat)
    end
    return _sum_of_square_norm*π/N^2
end

function volume_by_cross_section(model::AbstractCelegansModel; delta = 0.001)
    # integrate along the z-axis
    _central_spline = central_spline(model)
    last_z = _central_spline(0)[3]
    last_area = cross_section_area(model, 0)
    volume = 0
    for p in delta:delta:1
        z = _central_spline(p)[3]
        area = cross_section_area(model, p)
        volume += (area + last_area)/2 * (z-last_z)
        last_z = z
        last_area = area
    end
    return volume
end