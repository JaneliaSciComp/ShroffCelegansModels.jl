using GeometryBasics: Point3
using LinearAlgebra

# Single synthetic cross section at the origin: tangent = +z, right = +x.
function synthetic_contour(radii::AbstractVector)
    ShroffCelegansModels.FourierContour(
        [Point3(0.0, 0, 0)], [Point3(0.0, 0, 1)], [Point3(1.0, 0, 0)], reshape(float.(radii), :, 1))
end
inplane(θ, r = 1.0) = Point3(r * cos(θ), r * sin(θ), 0.0)

@testset "FourierContour synthetic cross sections" begin
    FC = ShroffCelegansModels
    # Circle: exact at every angle, true max equals the radius
    c = synthetic_contour(fill(7.0, 16))
    for θ in range(-π, π, length = 25)
        @test FC.contour_radius(c, 1, inplane(θ)) ≈ 7.0
    end
    @test FC.true_max_radius(c, 1) ≈ 7.0

    # Ellipse (semi-axes a=10 along x, b=6 along y) sampled at N angles
    a, b = 10.0, 6.0
    ell(θ) = a * b / sqrt((b * cos(θ))^2 + (a * sin(θ))^2)
    for N in (15, 16, 32)
        θs = 2π .* (0:N-1) ./ N
        c = synthetic_contour(ell.(θs))
        # exact at the sample angles
        for (k, θ) in enumerate(θs)
            @test FC.contour_radius(c, 1, inplane(θ)) ≈ ell(θ) atol = 1e-9
        end
        # close in between, and the max is the semi-major axis
        for θ in range(0, 2π, length = 37)
            @test FC.contour_radius(c, 1, inplane(θ)) ≈ ell(θ) rtol = 0.01
        end
        @test FC.true_max_radius(c, 1) ≈ a rtol = 0.005
        @test FC.true_max_radius(c, 1) >= maximum(ell.(θs)) - 1e-12
        # cached
        @test FC.true_max_radius(c, 1) == c.max_cache[1]
    end

    # Rotated ellipse (long axis at 45°): angle convention
    N = 16
    θs = 2π .* (0:N-1) ./ N
    c = synthetic_contour(ell.(θs .- π / 4))
    @test FC.contour_radius(c, 1, inplane(π / 4)) ≈ a atol = 1e-6
    @test FC.contour_radius(c, 1, inplane(π / 4 + π / 2)) ≈ b atol = 1e-6
end

@testset "FourierContour edge cases" begin
    FC = ShroffCelegansModels
    # odd/even N, small N
    for N in (1, 2, 3, 4, 5, 8)
        c = synthetic_contour(fill(3.0, N))
        @test FC.contour_radius(c, 1, inplane(0.3)) ≈ 3.0
        @test FC.true_max_radius(c, 1) ≈ 3.0
    end
    # even N result stays real/finite for uneven radii
    c = synthetic_contour([4.0, 2.0, 4.0, 2.0])
    @test isfinite(FC.contour_radius(c, 1, inplane(0.7)))
    @test FC.true_max_radius(c, 1) >= 4.0 - 1e-12

    c = synthetic_contour(fill(5.0, 8))
    # point on the central spline (zero vector)
    @test isfinite(FC.contour_radius(c, 1, Point3(0.0, 0, 0)))
    # point on the tangent axis (no in-plane component)
    @test isfinite(FC.contour_radius(c, 1, Point3(0.0, 0, 3.0)))
    # out-of-plane offset does not change the angle
    @test FC.contour_radius(c, 1, Point3(1.0, 1.0, 9.0)) ≈ 5.0
end

@testset "nearest_central_pt fourier_contour matches brute force" begin
    FC = ShroffCelegansModels
    model = FC.build_celegans_model(joinpath(@__DIR__, "fixtures", "lattice.csv"))
    lat = vec(FC.lattice(model))
    # lattice points plus deterministic offsets that probe the outer boundary
    offsets = [Point3(35.0, 0, 0), Point3(0.0, 25.0, 0), Point3(-20.0, 20.0, 10.0), Point3(60.0, -60.0, 0)]
    pts = vcat(lat, [p + o for p in lat[1:2:end] for o in offsets])
    efs = [1.0, 1.25, 1.5, 2.0, 3.0]

    z = LinRange(0, 1, length(model))
    cs = FC.central_spline(model)
    central_pts = cs.(z)
    c = FC.FourierContour(model, z)

    function reference(pt)
        ud = norm.(central_pts .- pt)
        for ef in efs
            d = copy(ud)
            for i in eachindex(z)
                d[i] <= FC.contour_radius(c, i, pt) * ef || (d[i] = Inf)
            end
            a = argmin(d)
            isfinite(d[a]) && return (z[a], d[a])
        end
        a = argmin(ud)
        return (z[a], Inf)
    end

    t, ncp, dist = FC.nearest_central_pt(model, pts, efs; fourier_contour = true)
    nfinite = 0
    for (k, pt) in enumerate(pts)
        zr, dr = reference(pt)
        @test t[k] == zr
        @test dist[k] == dr
        nfinite += isfinite(dr)
    end
    @test 0 < nfinite  # the comparison exercised the bounded branch
    # Pruning must be exact: the lazy true max can never undercut the sampled max
    mr = FC.max_radius_function(model).(z)
    @test all(FC.true_max_radius(c, i) >= mr[i] - 1e-9 for i in eachindex(z))

    # default path unchanged by the new keyword
    t0, _, d0 = FC.nearest_central_pt(model, pts, efs)
    t1, _, d1 = FC.nearest_central_pt(model, pts, efs; fourier_contour = false)
    @test t0 == t1 && d0 == d1
end

@testset "FourierContour orientation" begin
    FC = ShroffCelegansModels
    # Asymmetric profile: orient=-1 with radii sampled in the opposite sense
    # must reproduce the same physical contour.
    N = 16
    f(θ) = 6 + 2cos(θ) + 1.5sin(2θ) + 0.5sin(θ)
    θs = 2π .* (0:N-1) ./ N
    cpos = synthetic_contour(f.(θs))
    cneg = FC.FourierContour([Point3(0.0, 0, 0)], [Point3(0.0, 0, 1)], [Point3(1.0, 0, 0)],
        reshape(f.(-θs), :, 1), -1.0)
    for θ in range(-π, π, length = 21)
        @test FC.contour_radius(cpos, 1, inplane(θ)) ≈ f(θ) atol = 1e-6
        @test FC.contour_radius(cneg, 1, inplane(θ)) ≈ f(θ) atol = 1e-6
    end

    # Real model: orientation is ±1 and radial segment k lies near angle 2π(k-1)/N
    model = FC.build_celegans_model(joinpath(@__DIR__, "fixtures", "lattice.csv"))
    z = LinRange(0, 1, length(model))
    c = FC.FourierContour(model, z)
    @test abs(c.orient) == 1
    tss = FC.transverse_splines(model)
    Ns = length(tss)
    cs = FC.central_spline(model)
    dcs = FC.Derivative(1) * cs
    for i in (5, length(z) ÷ 2, length(z) - 5), k in 1:Ns
        v = tss[k](z[i]) - cs(z[i])
        θ = c.orient * FC._contour_angle(v, c.rights[i], c.tangents[i])
        expected = 2π * (k - 1) / Ns
        Δ = mod(θ - expected + π, 2π) - π
        @test abs(Δ) < 0.3
    end
end

@testset "FourierContour true_max_radius robustness and cache" begin
    FC = ShroffCelegansModels
    # deterministic pseudo-random radii, including high N and spiky profiles
    for N in (4, 7, 16, 31, 32), seed in 1:6
        radii = [5 + 3 * sin(k^2 * 1.7 + seed * 2.3) + 1.2 * cos(k * seed) for k in 1:N]
        c = synthetic_contour(radii)
        M = 20000
        brute = maximum(FC._contour_eval(c, 1, 2π * j / M) for j in 0:M-1)
        @test FC.true_max_radius(c, 1) ≈ brute rtol = 1e-6
        @test FC.true_max_radius(c, 1) >= maximum(radii) - 1e-9
    end

    # cache: starts empty, filled only for queried index, stable on repeat
    z = [0.0, 0.5, 1.0]
    c = FC.FourierContour(fill(Point3(0.0, 0, 0), 3), fill(Point3(0.0, 0, 1), 3),
        fill(Point3(1.0, 0, 0), 3), [5.0 6.0 7.0; 5.5 6.5 7.5; 5.0 6.0 7.0; 4.5 5.5 6.5])
    @test all(isnan, c.max_cache)
    v = FC.true_max_radius(c, 2)
    @test isfinite(v) && isnan(c.max_cache[1]) && isnan(c.max_cache[3])
    @test c.max_cache[2] == v
    @test FC.true_max_radius(c, 2) === v
end

@testset "nearest_central_pt consistency across expansion factors" begin
    FC = ShroffCelegansModels
    model = FC.build_celegans_model(joinpath(@__DIR__, "fixtures", "lattice.csv"))
    lat = vec(FC.lattice(model))
    pts = vcat(lat, [p + Point3(30.0, 10.0, 0) for p in lat[1:3:end]], [Point3(1e4, 1e4, 1e4)])
    efs = [1.0, 1.5, 2.0, 3.0]
    ds = [last(FC.nearest_central_pt(model, pts, ef; fourier_contour = true)) for ef in efs]
    for j in 1:length(efs)-1, k in eachindex(pts)
        # a larger expansion factor admits a superset of candidates
        @test ds[j+1][k] <= ds[j][k]
    end
    # far point finds nothing and does not throw
    @test all(isinf, getindex.(ds, length(pts)))
    # a sequence of factors yields the first factor with any finite candidate
    _, _, dseq = FC.nearest_central_pt(model, pts, efs; fourier_contour = true)
    for k in eachindex(pts)
        j = findfirst(d -> isfinite(d[k]), ds)
        @test dseq[k] == (j === nothing ? Inf : ds[j][k])
    end
end

@testset "default path matches original algorithm" begin
    FC = ShroffCelegansModels
    model = FC.build_celegans_model(joinpath(@__DIR__, "fixtures", "lattice.csv"))
    lat = vec(FC.lattice(model))
    pts = vcat(lat, [p + Point3(30.0, 10.0, 0) for p in lat[1:3:end]])
    efs = [1.0, 1.25, 1.5, 2.0, 3.0]
    z = LinRange(0, 1, length(model))
    central_pts = FC.central_spline(model).(z)
    radius = FC.max_radius_function(model).(z)
    # verbatim copy of the pre-change algorithm
    function original(pt)
        ud = norm.(central_pts .- pt)
        a = 1
        dist = copy(ud)
        for ef in efs
            dist[ud .> radius .* ef] .= Inf
            a = argmin(dist)
            dist[a] != Inf && break
            dist = copy(ud)
        end
        z[a], central_pts[a], dist[a]
    end
    t, ncp, d = FC.nearest_central_pt(model, pts, efs)
    for (k, pt) in enumerate(pts)
        zo, co, do_ = original(pt)
        @test t[k] == zo
        @test ncp[k] == co
        @test d[k] == do_
    end
    # untwist_annotations is unchanged by the keyword's presence
    @test FC.untwist_annotations(model, pts) == FC.untwist_annotations(model, pts; fourier_contour = false)
end

@testset "prebuilt FourierContour reuse" begin
    FC = ShroffCelegansModels
    model = FC.build_celegans_model(joinpath(@__DIR__, "fixtures", "lattice.csv"))
    lat = vec(FC.lattice(model))
    pts = vcat(lat, [p + Point3(30.0, 10.0, 0) for p in lat[1:3:end]])
    efs = [1.0, 1.5, 2.0, 3.0]
    z = LinRange(0, 1, length(model))
    c = FC.FourierContour(model, z)
    fresh = FC.nearest_central_pt(model, pts, efs; fourier_contour = true)
    reused = FC.nearest_central_pt(model, pts, efs; fourier_contour = c)
    @test fresh == reused
    # cache is populated by the first use and reused by later calls
    @test any(!isnan, c.max_cache)
    snapshot = copy(c.max_cache)
    again = FC.nearest_central_pt(model, pts[1:5], efs; fourier_contour = c)
    # previously cached maxima are retained unchanged
    @test all(c.max_cache[i] == snapshot[i] for i in eachindex(snapshot) if !isnan(snapshot[i]))
    @test again == FC.nearest_central_pt(model, pts[1:5], efs; fourier_contour = true)
    # untwist paths accept it too
    @test FC.untwist_annotations(model, pts; fourier_contour = c) ==
          FC.untwist_annotations(model, pts; fourier_contour = true)
    @test FC.untwist_annotation(model, pts[1]; fourier_contour = c) ==
          FC.untwist_annotation(model, pts[1]; fourier_contour = true)
    # mismatched contour is rejected
    bad = FC.FourierContour(fill(Point3(0.0, 0, 0), 3), fill(Point3(0.0, 0, 1), 3),
        fill(Point3(1.0, 0, 0), 3), fill(5.0, 8, 3))
    @test_throws ArgumentError FC.nearest_central_pt(model, pts, efs; fourier_contour = bad)
end
