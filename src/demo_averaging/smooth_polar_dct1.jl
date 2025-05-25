using CoordinateTransformations
using FFTW: fftfreq, fft, ifft

function smooth_polar_dct1(positions_over_time, σ_r, σ_θ, σ_z = 0)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()

    N = length(positions_over_time)

    # DCT Type I Mirroring
    # positions_over_time = [positions_over_time; positions_over_time[end-1:-1:2]]

    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _r = [_r; _r[end-1:-1:2]]
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    if σ_r > 0
        _filter = G(σ_r, length(_r))
        _r = abs.(ifft(fft(_r) .* _filter))
    end
    _r = @view _r[1:N]

    _cs = exp.(_θ .* 1im)
    _cs = [_cs; _cs[end-1:-1:2]]
    if σ_θ > 0
        _filter = G(σ_θ, length(_cs))
        _cs = ifft(fft(_cs) .* _filter)
    end
    _θ = angle.(_cs)

    _z = [_z; _z[end-1:-1:2]]
    if σ_z > 0
        lpz_filter = G(0.5, length(_z))
        lp_z = real.(ifft(fft(_z) .* lpz_filter))
        _filter = G(σ_z, length(_z))
        _z = real.(ifft(fft(_z .- lp_z) .* _filter))
        _z .+= lp_z
    end
    _z = @view _z[1:N]

    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
    # return @view positions_over_time[1:N]
end