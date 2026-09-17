using ShroffCelegansModels.CoordinateTransformations
using ShroffCelegansModels.FFTW: fftfreq, fft, ifft

_replicate_pad(x, pad::Int) = pad > 0 ? [fill(x[1], pad); x; fill(x[end], pad)] : x

function smooth_polar_dct1(positions_over_time, σ_r, σ_θ, σ_z = 0; edge_pad::Int = 0)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()

    N = length(positions_over_time)
    lo, hi = edge_pad + 1, edge_pad + N

    # DCT Type I Mirroring
    # positions_over_time = [positions_over_time; positions_over_time[end-1:-1:2]]

    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _r = _replicate_pad(_r, edge_pad)
    _r = [_r; _r[end-1:-1:2]]
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    if σ_r > 0
        _filter = G(σ_r, length(_r))
        _r = abs.(ifft(fft(_r) .* _filter))
    end
    _r = @view _r[lo:hi]

    _cs = exp.(_θ .* 1im)
    _cs = _replicate_pad(_cs, edge_pad)
    _cs = [_cs; _cs[end-1:-1:2]]
    if σ_θ > 0
        _filter = G(σ_θ, length(_cs))
        _cs = ifft(fft(_cs) .* _filter)
    end
    _θ = @view angle.(_cs)[lo:hi]

    _z = _replicate_pad(_z, edge_pad)
    _z = [_z; _z[end-1:-1:2]]
    if σ_z > 0
        lpz_filter = G(0.5, length(_z))
        lp_z = real.(ifft(fft(_z) .* lpz_filter))
        _filter = G(σ_z, length(_z))
        _z = real.(ifft(fft(_z .- lp_z) .* _filter))
        _z .+= lp_z
    end
    _z = @view _z[lo:hi]

    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
end
