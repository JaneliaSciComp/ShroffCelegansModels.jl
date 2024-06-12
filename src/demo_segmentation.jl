using FileIO
using TiffImages

const seg_image_path = raw"\\nearline4.hhmi.org\shroff\shrofflab\NK2436_pat-3_mNG\NK2436_Pat-3"
const seg_image_path_pos0_spimb = joinpath(seg_image_path, "Pos0", "SPIMB")

function load_img_dir(
    dirpath::AbstractString = pwd(),
    timepoints = nothing)
    if isnothing(timepoints)
        timepoints = extract_timepoints_from_files(readdir(dirpath))
    end
    map(timepoints) do tp
        filename = "SPIMB-$tp.tif"
        img = load(filename)[:,1:1024,:]
        maximum(img; dims=3)
    end
end
function extract_timepoints_from_files(filenames::Vector{<: AbstractString})
    timepoints = map(filenames) do filename
        m = match(r"([0-9]+).tif$", filename)
        num = tryparse(Int, m.captures[1])
        isnothing(num) ? -1 : num
    end
    sort(filter(>=(0), timepoints))
end

function gaussian2d(sigma=3, extent=sigma * 3)
    F = map(CartesianIndices((-extent:extent, -extent:extent))) do ci
        t = ci.I
        exp(-sum(t .^ 2) / 2 / sigma^2)
    end
    return F ./ sum(F)
end
function gaussian2d_deriv1(sigma=3, extent=sigma * 3; dim=1)
    sigma2 = sigma^2
    F = map(CartesianIndices((-extent:extent, -extent:extent))) do ci
        t = ci.I
        -t[dim] / sigma2 * exp(-sum(t .^ 2) / 2 / sigma2)
    end
end
function gaussian2d_deriv2(sigma=3, extent=sigma * 3; dim=1)
    sigma2 = sigma^2
    sigma4 = sigma^4
    F = map(CartesianIndices((-extent:extent, -extent:extent))) do ci
        t = ci.I
        -(t[dim]^2 - sigma2)/sigma4 * exp(-sum(t .^ 2) / 2 / sigma2)
    end
end
function gaussian2d_deriv2_diag(sigma=3, extent=sigma * 3)
    sigma2 = sigma^2
    sigma4 = sigma^4
    F = map(CartesianIndices((-extent:extent, -extent:extent))) do ci
        t = ci.I
        (t[1]*t[2])/sigma4 * exp(-sum(t .^ 2) / 2 / sigma2)
    end
end

