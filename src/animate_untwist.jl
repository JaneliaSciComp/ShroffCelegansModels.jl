using CoordinateTransformations
using GeometryBasics
using LinearAlgebra

function stretch_worm(pts, faces, spts, e)
    pts0 = pts .- pts[2]
    central_pts = @view pts0[2:3:end]
    spts = ((1 - e) .* central_pts) .+ (e .* spts)
    pts0[1:3:end] .-= central_pts
    pts0[1:3:end] .+= spts
    pts0[3:3:end] .-= central_pts
    pts0[3:3:end] .+= spts
    central_pts .= spts
    return GeometryBasics.Mesh(swapyz.(pts0), faces), swapyz.(pts0)
end

function untwist_worm(pts, faces, spts, e)
    pts0 = pts .- pts[2]

    L_pts = @view pts0[1:3:end]
    central_pts = @view pts0[2:3:end]
    R_pts = @view pts0[3:3:end]

    L_pts .-= central_pts
    R_pts .-= central_pts

    L_norm = norm.(L_pts)
    R_norm = norm.(R_pts)

    L_unit_vec = (1-e) .* (L_pts ./ L_norm) .+ (e .* Point3f( 1, 0, 0),)
    R_unit_vec = (1-e) .* (R_pts ./ R_norm) .+ (e .* Point3f(-1, 0, 0),)

    sfc = SphericalFromCartesian()
    cfs = CartesianFromSpherical()
    
    L_sph = sfc.(L_unit_vec)
    R_sph = sfc.(R_unit_vec)

    #L_pts .= cfs.(Spherical.(L_norm, L_sph.θ, L_sph.ϕ))
    #R_pts .= cfs.(Spherical.(R_norm, R_sph.θ, L_sph.ϕ))

    L_pts .= map(L_norm, L_sph) do L_norm, L_sph
        cfs(Spherical(L_norm, L_sph.θ, L_sph.ϕ))
    end

    R_pts .= map(R_norm, R_sph) do R_norm, R_sph
        cfs(Spherical(R_norm, R_sph.θ, R_sph.ϕ))
    end

    #=
    L_pts .*= (1-e)
    R_pts .*= (1-e)

    L_pts .+= e .* Point3.(L_norm, 0.0, 0.0)
    R_pts .-= e .* Point3.(R_norm, 0.0, 0.0)
    =#

    L_pts .+= spts
    central_pts .= spts
    R_pts .+= spts

    return GeometryBasics.Mesh(swapyz.(pts0), faces), swapyz.(pts0)
    second = x -> x[2]
    theta = atan.(second.(pts0), first.(pts0))
    gamma = atan.(last.(pts0), first.(pts0))
    return theta, L_norm, R_norm, pts0
    L_theta = @view(theta[1:3:end])
    R_theta = @view(theta[3:3:end])
    L_theta = 0
    R_theta = π
    z = last.(spts)
    pts0[1:3:end] .= Point3.(cos.(L_theta) .* L_norm, sin.(L_theta) .* L_norm, z)
    pts0[3:3:end] .= Point3.(cos.(R_theta) .* R_norm, sin.(R_theta) .* R_norm, z)
    #return GeometryBasics.Mesh(swapyz.(pts0), faces)
    return theta
end

function untwist_worm_spherical(pts, faces, spts, e, f)
    pts0 = pts .- pts[1]
    L_pts = @view pts0[1:3:end]
    central_pts = @view pts0[2:3:end]
    R_pts = @view pts0[3:3:end]
    L_pts .-= central_pts
    R_pts .-= central_pts
    sfc = SphericalFromCartesian()
    cfs = CartesianFromSpherical()
    L_pts .= map(sfc.(L_pts)) do sph_pt
        cfs(Spherical(sph_pt.r, sph_pt.θ * (1-e), sph_pt.ϕ *(1-f)))
    end
    R_pts .= map(sfc.(R_pts)) do sph_pt
        cfs(Spherical(sph_pt.r, sph_pt.θ * (1-e) + e * π, sph_pt.ϕ *(1-f)))
    end
    pts0[1:3:end] .+= spts
    pts0[3:3:end] .+= spts
    central_pts .= spts
    return GeometryBasics.Mesh(swapyz.(pts0), faces)
    second = x -> x[2]
    theta = atan.(second.(pts0), first.(pts0))
    gamma = atan.(last.(pts0), first.(pts0))
    return theta, L_norm, R_norm, pts0
    L_theta = @view(theta[1:3:end])
    R_theta = @view(theta[3:3:end])
    L_theta = 0
    R_theta = π
    z = last.(spts)
    pts0[1:3:end] .= Point3.(cos.(L_theta) .* L_norm, sin.(L_theta) .* L_norm, z)
    pts0[3:3:end] .= Point3.(cos.(R_theta) .* R_norm, sin.(R_theta) .* R_norm, z)
    #return GeometryBasics.Mesh(swapyz.(pts0), faces)
    return theta
end

function animate_untwist(smodel; fig = nothing)
    model = smodel.twisted_model
    _pts, _faces = ShroffCelegansModels.get_model_manifold_mesh_components(model)
    r = LinRange(0,1,length(smodel))
    spts = ShroffCelegansModels.central_spline(smodel).(r)

    # fig = Figure()

    # ax = Axis3(fig[1,1], aspect=:data)

    cs_pts = ShroffCelegansModels.central_spline(model).(r)
    scs_pts = ShroffCelegansModels.central_spline(smodel).(r)

    if isnothing(fig)
        out = lines([cs_pts; swapyz.(scs_pts)]; visible = false)
        fig = out.figure
        ax = out.axis
        resize!(fig, 1920, 1080)
        display(fig)
    else
        ax = fig.current_axis[]
    end
    #return fig

    #display(fig)
    
    M, pts = stretch_worm(_pts, _faces, spts, 0.0)
    M = Observable(M)
    L = Observable(@view(pts[1:3:end]))
    C = Observable(@view(pts[2:3:end]))
    R = Observable(@view(pts[3:3:end]))

    model_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points = p->swapyz(p - cs_pts[1])))
    mesh_plot = mesh!(ax, model_mesh; transparency = true)
    display(fig)
    sleep(1)

    for rep in 1:3
        mesh!(ax, M; shading = MakieCore.automatic, color = [i for c in 1:length(model) for i in 1:3], colormap = :buda)
        sleep(1)

        for i in 0:10:100
            mesh_plot.alpha[] = (100-i)/100
            sleep(0.1)
        end

        #return fig
        lines!(ax, L, color=:red)
        lines!(ax, C, color=:magenta)
        lines!(ax, R, color=:green)

        display(fig)

        for i=0:10:100
            M[], pts = stretch_worm(_pts, _faces, spts, 0.01*i)
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
            if i == 0
                sleep(1)
            else
                sleep(0.1)
            end
        end
        for i=0:100
            M[], pts = untwist_worm(_pts, _faces, spts, 0.01*i)
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
            if i == 0
                sleep(0.5)
            else
                sleep(0.01)
            end
        end
        sleep(2)
        for i=100:-1:0
            M[], pts = untwist_worm(_pts, _faces, spts, 0.01*i)
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
            if i == 0
                sleep(0.5)
            else
                sleep(0.01)
            end
        end
        for i=100:-10:0
            M[], pts = stretch_worm(_pts, _faces, spts, 0.01*i)
            L[] = @view(pts[1:3:end])
            C[] = @view(pts[2:3:end])
            R[] = @view(pts[3:3:end])
            if i == 0
                sleep(1)
            else
                sleep(0.1)
            end
        end

        for i in 100:-10:0
            mesh_plot.alpha[] = (100-i)/100
            sleep(0.1)
        end

    end
    return fig
end