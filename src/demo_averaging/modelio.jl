using HDF5
using Printf
using BSplineKit: SplineInterpolation
using GeometryBasics: Point3f, Point3
using GLMakie

function save_celegans_model(
    parent::Union{HDF5.File, HDF5.Group},
    name::String,
    model::ShroffCelegansModels.CelegansModel
)
    g = create_group(parent, name)
    # g["names"] = model.names
    h5_string3 = HDF5.datatype(eltype(model.names))
    names_ds = create_dataset(g, "names", h5_string3, (length(model.names),))
    write_dataset(names_ds, h5_string3, model.names)
    ##names_ds[1] = model.names[1]
    # central_spline_g = create_group(g, "central_spline")
    central_spline_g = save_spline_interpolation(
        g,
        "central_spline",
        model.central_spline,
        "Central B-spline that goes through the midpoints between seam cells of the C. elegans model."
    )
    # central_spline_g["order"] = BSplineKit.order(model.central_spline)
    transverse_splines_g = create_group(g, "transverse_splines")
    for (i, s) in pairs(model.transverse_splines)
        save_spline_interpolation(
            transverse_splines_g,
            @sprintf("transverse_spline_%02d",i),
            s,
            "Transverse B-spline along the contour from anterior to posterior of the C. elegans model"
        )
    end
    a = attrs(g)
    a["type_description"] = "C. elegans model"
    a["julia_type"] = string(typeof(model))
end

function save_spline_interpolation(
    parent::Union{HDF5.File, HDF5.Group},
    name::String,
    spline::SplineInterpolation,
    description::String = ""
)
    pts = BSplineKit.SplineInterpolations.interpolation_points(spline)
    collocation_points = spline.(pts)
    order = BSplineKit.order(spline)
    g = create_group(parent, name)
    g["abscissa"] = pts
    g["ordinate"] = collocation_points
    g["coefficients"] = coefficients(spline)
    a = attrs(g)
    a["order"] = order
    a["boundary_conditions"] = "Natural"
    order_name =
        order == 1 ? "constant" :
        order == 2 ? "linear" :
        order == 3 ? "quadratic" :
        order == 4 ? "cubic" :
        string(order)
    # TODO: Detect boundary conditions
    a["type_description"] = "B-Spline interpolation of order $(order) ($(order_name)) with natural boundary conditions"
    a["description"] = description
    a["julia_type"] = string(typeof(spline))
    return g
end

function load_spline_interpolation(g::Union{HDF5.File, HDF5.Group})
    x = read(g["abscissa"])
    y = read(g["ordinate"], Point3{Float64})
    a = attrs(g)
    a["order"] == 4 || error("Only order 4 (Cubic) is implemented for loading.")
    a["boundary_conditions"] == "Natural" || error("Only Natural boundary conditions are implemented")
    return ShroffCelegansModels.interpolate_natural_cubic_spline(x, y)
end

function load_celegans_model(g::Union{HDF5.File, HDF5.Group})
    central_spline = load_spline_interpolation(g["central_spline"])
    transverse_splines = map(1:32) do i
        load_spline_interpolation(g[@sprintf("transverse_splines/transverse_spline_%02d", i)])
    end
    names = g["names"][]
    return ShroffCelegansModels.Types.CelegansModel(transverse_splines, central_spline, names)
end

