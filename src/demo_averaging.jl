using BSplineKit.SplineInterpolations: interpolation_points
using ShroffCelegansModels
using ShroffCelegansModels.Types: AbstractCelegansModel
using JSON3
using Missings
using StatsBase
includet("makie.jl")
using Makie
using CSV
using DataFrames
using BSplineKit
using ThinPlateSplines
using CoordinateTransformations
using FFTW
using ProgressMeter

include("demo_averaging/read_config_json.jl")
include("demo_averaging/get_lattice.jl")
include("demo_averaging/build_models_over_time.jl")
include("demo_averaging/average_sliders.jl")
include("demo_averaging/cross_sections_at_knots.jl")
include("demo_averaging/seam_cell_pts.jl")
include("demo_averaging/get_avg_models.jl")
include("demo_averaging/show_average_models.jl")


include("demo_averaging/get_straightened_annotations.jl")
include("demo_averaging/load_straightened_annotations_over_time.jl")

using ShroffCelegansModels.MIPAVIO: get_model_csv, get_integrated_annotations,
    get_straightened_lattice, get_straightened_lattice_xy_center,
    mipav_df_to_points, mipav_df_to_points_dict
include("demo_averaging/get_cell_trajectory_dict.jl")
include("demo_averaging/transform_annotations.jl")
include("demo_averaging/plot_int_cells.jl")
include("demo_averaging/show_average_models_with_annotations.jl")
include("demo_averaging/show_average_models_with_annotations_demo.jl")
include("demo_averaging/NormalizedTimeFunction.jl")
include("demo_averaging/debug_average_models_with_annotations.jl")
include("demo_averaging/show_average_annotations.jl")
include("demo_averaging/check_dataset.jl")
include("demo_averaging/modelio.jl")

#=
function HDF5.hdf5_type_id(::Type{String3})
    HDF5.API.h5t_create(HDF5.API.H5T_STRING,4);
end
function HDF5.datatype(t::Type{String3})
    return HDF5.Datatype(HDF5.hdf5_type_id(t))
end
=#

include("demo_averaging/check_annotation_within_contour.jl")
# include("demo_averaging/loading.jl")
include("demo_averaging/stretched_analysis.jl")
include("demo_averaging/angle_heatmap.jl")
