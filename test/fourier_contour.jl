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
