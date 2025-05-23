using ShroffCelegansModels: AbstractCelegansModel, transverse_splines, get_model_contour_mesh, num_transverse_splines
using ColorSchemes: colorschemes
using MakieCore: MakieCore
using GeometryBasics: Point3, Point3f
using ShroffCelegansModels.Points: cross_sections
using Observables: throttle

function plot_celegans_model(model::AbstractCelegansModel)
    r = LinRange(0, 1, length(model))
    splines = transverse_splines(model)
    h = lines(splines[1].(r))
    for i in 2:32
        lines!(splines[i].(r))
    end
    return h
end

function plot_celegans_model!(model::AbstractCelegansModel)
    r = LinRange(0, 1, length(model))
    splines = transverse_splines(model)
    for i in 1:32
        lines!(splines[i].(r))
    end
end

const cyclic_colorschemes = [
    :cyclic_grey_15_85_c0_n256,
    :cyclic_grey_15_85_c0_n256_s25,
    :cyclic_mrybm_35_75_c68_n256,
    :cyclic_mrybm_35_75_c68_n256_s25,
    :cyclic_mygbm_30_95_c78_n256,
    :cyclic_mygbm_30_95_c78_n256_s25,
    :cyclic_protanopic_deuteranopic_bwyk_16_96_c31_n256,
    :cyclic_protanopic_deuteranopic_wywb_55_96_c33_n256,
    :cyclic_tritanopic_cwrk_40_100_c20_n256,
    :cyclic_tritanopic_wrwc_70_100_c20_n256,
    :cyclic_wrwbw_40_90_c42_n256,
    :cyclic_wrwbw_40_90_c42_n256_s25,
]


function MakieCore.mesh(
    model::AbstractCelegansModel;
    colorscheme = :cyclic_tritanopic_cwrk_40_100_c20_n256,
    shading = MakieCore.automatic,
    n_ellipse_pts = num_transverse_splines(model),
    color=repeat(colorschemes[colorscheme][1:256÷n_ellipse_pts:256], length(model)),
    colorrange = (1,n_ellipse_pts),
    alpha = 0.5,
    transform_points = identity,
    kwargs...
)
    worm_mesh = get_model_contour_mesh(model; transform_points)
    mesh(
        worm_mesh;
        shading,
        color,
        colorrange,
        alpha,
        kwargs...
    )
end

function MakieCore.mesh!(
    axis,
    model::AbstractCelegansModel;
    colorscheme = :cyclic_wrwbw_40_90_c42_n256,
    shading = MakieCore.automatic,
    n_ellipse_pts = num_transverse_splines(model),
    color=repeat(colorschemes[colorscheme][1:256÷n_ellipse_pts:256], length(model)),
    colorrange = (1,n_ellipse_pts),
    alpha = 0.7,
    transform_points = identity,
    kwargs...
)
    worm_mesh = get_model_contour_mesh(model; transform_points)
    mesh!(
        axis,
        worm_mesh;
        shading,
        color,
        colorrange,
        alpha,
        kwargs...
    )
end

MakieCore.convert_arguments(P::Type{<: MakieCore.Mesh}, model::AbstractCelegansModel) =
    MakieCore.convert_arguments(P, get_model_contour_mesh(model))
MakieCore.convert_arguments(P::Type{<: MakieCore.Mesh}, model::AbstractCelegansModel, f::Function) =
    MakieCore.convert_arguments(P, get_model_contour_mesh(model, transform_points = f))

function MakieCore.convert_arguments(P::Type{<: MakieCore.Lines}, model::AbstractCelegansModel, f::Function = identity)
    splines = transverse_splines(model)
    r = LinRange(0, 1, length(model))
    pts = map(splines) do spline
        spline.(r)
    end
    #=
    gaps = fill!(similar(pts[1]), Point3f(NaN))
    cross_sections = hcat(pts..., pts[1], gaps)
    cross_sections = permutedims(cross_sections, (2,1)) |> vec
    =#

    gaps = fill!(similar(pts), [Point3(NaN,NaN,NaN)])
    pts = permutedims(hcat(pts, gaps), [2,1]) |> vec
    pts = vcat(pts...)
    pts = f.(pts)
    MakieCore.convert_arguments(P, pts)
end

#function MakieCore.convert_arguments(P::Type{<: MakieCore.Text}, model::AbstractCelegansModel, f::Function = identity)
function MakieCore.text!(ax, model::AbstractCelegansModel, f::Function = identity)
    pts = ShroffCelegansModels.interpolation_points(model)
    splines = transverse_splines(model)
    # This probably swapped. 1 should be right, 17 should be left
    left_spline = splines[1]
    right_spline = splines[17]
    text!(ax, f.(left_spline.(pts)); text = ShroffCelegansModels.Types.names(model)[1:2:end])
    # text!(ax, f.(right_spline.(pts)); text = ShroffCelegansModels.Types.names(model)[2:2:end])
end

# const voxel_size = 0.1625 # um

swapyz(p::P) where P <: Point3 = P(p[1], p[3], p[2])
swapyz_scale(p::P) where P <: Point3 = P(p[1], p[3], p[2]) * voxel_size
swapyz_unscale(p::P) where P <: Point3 = P(p[1], p[3], p[2]) / voxel_size
scale_dims(p::P) where P <: Point3 = P(p[1], p[2], p[3]) * voxel_size
