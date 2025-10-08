module ShroffCelegansModels
    using BSplineKit
    using CSV
    using ColorSchemes: ColorSchemes
    using Colors: Colors
    using CoordinateTransformations: CoordinateTransformations
    using DataFrames
    using Dates
    using FFTW
    using FileIO: FileIO
    using FixedPointNumbers: FixedPointNumbers
    using GeometryBasics
    using HDF5: HDF5
    using Interpolations: Interpolations
    using JSON3: JSON3
    using LRUCache
    using LinearAlgebra
    using Makie
    using Missings: Missings
    using Observables: Observables
    using Pkg: Pkg
    using PrecompileTools: @setup_workload, @compile_workload
    using ProgressMeter: ProgressMeter
    using QuadGK
    using Sockets: Sockets
    using StaticArrays: StaticArrays
    using Statistics: Statistics
    using StatsBase: StatsBase
    using ThinPlateSplines: ThinPlateSplines
    using TiffImages: TiffImages
    using ShroffCelegansModelsCore:
        ShroffCelegansModelsCore,
        get_circle_points,
        build_celegans_model,
        get_model_contour_mesh,
        get_sections,
        get_model_manifold_mesh_components,
        straighten_celegans_model,
        psuedo_seam_cells,
        average,
        upsample,
        ModelTimeSeries,
        StraightenedModelTimeSeries,
        get_lattice,
        get_lattice_filepath,
        nextModelIndex,
        prevModelIndex,
        radial_cross_section,
        cross_section_area,
        volume_by_cross_section,
        nearest_central_pt,
        get_central_point_parameters,
        max_radius_function,
        nearest_central_plane,
        untwist_annotations,
        untwist_annotation,
        twisted_annotations,
        distance_to_twisted_annotation

    using ShroffCelegansModelsCore.Datasets
    using ShroffCelegansModelsCore.MIPAVIO: MIPAVIO
    using ShroffCelegansModelsCore.ParametricSplines
    using ShroffCelegansModelsCore.Types
    using ShroffCelegansModelsCore.Points

    # Straightened annotations
    const annotations_cache = Dict{Tuple{String, UnitRange, Bool}, Vector}()
    # Warped annotations, with MIPAV straightening
    const annotation_position_cache = Dict{String, Any}()
    # Warped annotations, with Mark's straightening
    const my_annotation_position_cache = Dict{String, Vector{Vector{Point3{Float64}}}}()



    include("parse_worm_dataset_path.jl")

end
