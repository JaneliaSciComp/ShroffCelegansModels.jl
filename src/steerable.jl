using CoordinateTransformations
using StaticArrays
using GeometryBasics

function getFrequencySpaceCoordinates(Nx, Ny = Nx)
    hy = Ny÷2
    hx = Nx÷2
    ry = (0:Ny-1) .- hy
    rx = (0:Nx-1) .- hx
    cis = CartesianIndices((ry, rx))
    pfc = PolarFromCartesian()
    map(cis) do ci
        pt = Point2(Tuple(ci) ./ (hy, hx))
        pfc(pt)
    end
end

function angularKernel(K, angle = nothing, N = 1024)

    n = 2K + 1
    if isnothing(angle)
        # angle = (0:n-1)*(pi/n)
        angle = 0
    end

    coords = getFrequencySpaceCoordinates(N)
    
    θ = map(coords) do coord
        θ = coord.θ - angle
        mod(θ + π, 2π) - π
    end

    s_a = π/(2K+1)

    θ_s = θ ./ s_a

    angularFilter = 2*exp.(-θ_s.^2/2)
    angularFilter_shifted = @view angularFilter[[1; end:-1:2], [1; end:-1:2]]
    filterKernel = 0.5*(angularFilter .+ angularFilter_shifted)
    posMask = abs.(θ) .< π/2
    filterKernel = filterKernel .* (1 .+ 1im.*(posMask .*2 .- 1))

    return filterKernel
end
angularKernel(; K = 5, angle = nothing, N = 1024) = angularKernel(K, angle, N)

function radialKernel(f_c, b_f = nothing, N = 1024)
    if isnothing(b_f) 
        b_f = f_c / sqrt(2)
    end

    coords = getFrequencySpaceCoordinates(N)
    f = map(coords) do coord
        coord.r
    end

    if f_c == 0
        radialFilter = f.^2 ./ (2*b_f.^2)
        radialFilter[1] = 0.0
    else
        K_f = (f_c ./ b_f).^2
        f_s = f ./ f_c
        radialFilter = f_s.^K_f
        radialFilter = radialFilter .* exp.((1 .- f_s.^2) .* (K_f ./2))
    end
end
radialKernel(; f_c, b_f = nothing, N = 1024) = radialKernel(f_c, b_f, N)

using Interpolations
function nonMaximumSupression!(res, th, sup = zero(eltype(res)))
    Base.require_one_based_indexing(res)
    ny, nx = size(res)
    itp = interpolate((1:ny, 1:nx), res, Gridded(Linear()))
    etp = extrapolate(itp, Reflect())

    if size(th) == ()
        offset = reverse(sincos(th))
        map(CartesianIndices(res)) do ci
            A1 = etp((Tuple(ci) .+ offset)...)
            A2 = etp((Tuple(ci) .- offset)...)
            r = res[ci]
            if r < A1 || r < A2
                res[ci] = sup
            end
        end
    elseif size(th) == size(res)
        map(CartesianIndices(res), th) do ci, θ
            offset = reverse(sincos(θ))
            A1 = etp((Tuple(ci) .+ offset)...)
            A2 = etp((Tuple(ci) .- offset)...)
            r = res[ci]
            if r < A1 || r < A2
                res[ci] = sup
            end
        end
    else
        throw(ArgumentError("The size of `res`, $(size(res)), and `th`, $(size(th)), do not match."))
    end
    return res
end
function nonMaximumSupression(res, th, sup = zero(eltype(res)))
    return nonMaximumSupression!(copy(res), th, sup)
end
function nonLocalMaximaSupression(
    rotationResponse,
    theta = nothing,
    suppresionValue = zero(eltype(rotationResponse))
)
    nO = size(rotationResponse, 3)

    if isnothing(theta)
        theta = 0:nO - 1
        theta = theta .* π ./ nO
    end

    nlms = copy(rotationResponse)
    for o in 1:nO
        nlms[:,:,o] .= nonMaximumSupression(rotationResponse[:,:,o], theta[o], suppresionValue)
    end
    map(CartesianIndices(nlms)) do ci
        prev = mod1(ci[3]-1, nO)
        prev = CartesianIndex(ci[1], ci[2], prev)
        next = mod1(ci[3]+1, nO)
        next = CartesianIndex(ci[1], ci[2], next)
        r = rotationResponse[ci]
        if r < rotationResponse[next] ||
           r < rotationResponse[prev]
            suppresionValue
        else
            nlms[ci]
        end
    end
end