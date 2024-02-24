using BSplineKit.SplineInterpolations: interpolation_points
using ShroffCelegansModels
using ShroffCelegansModels.Types: AbstractCelegansModel
using JSON3
using Missings
using StatsBase
includet("makie.jl")
using GLMakie
using CSV
using DataFrames
using BSplineKit
using ThinPlateSplines

#const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_full.json"
const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_2024_01_10_v2.json"
const voxel_size = 0.1625 # um

function read_config_json(config_path::AbstractString = config_path)
    config_json = JSON3.read(config_path)
    cell_keys = Dict{String, Vector{ShroffCelegansModels.CellKey}}()
    map(config_json.data.strains) do strain
        cell_keys[strain.name] = map(strain.folderpaths) do folder_path
            cell_key_path = joinpath(folder_path, "cell_key.json")
            try
                cell_key_json = JSON3.read(cell_key_path)
                ShroffCelegansModels.CellKey(cell_key_json)
            catch err
                throw(ArgumentError("Problem parsing $cell_key_path"))
            end
        end
    end

    datasets = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
    map(config_json.data.strains) do strain
        datasets[strain.name] = map(strain.folderpaths) do folder_path
            ShroffCelegansModels.NormalizedDataset(joinpath(folder_path, "RegB"))
        end
    end
    return config_json, cell_keys, datasets
end

config_json, cell_keys, datasets = read_config_json()

flattened_datasets = collect(Iterators.flatten(values(datasets)))

min_length = map(flattened_datasets) do ds
    length(range(ds.cell_key))
end |> minimum

function get_lattice(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)::Union{Missing, String}
    timepoint = range(ds.cell_key)[time_offset]
    if timepoint ∈ ds.cell_key.outliers
        return missing
    else
        filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
        if isfile(filepath)
            return filepath
        else
            throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
        end
    end
end

function build_models_over_time(dataset::ShroffCelegansModels.Datasets.NormalizedDataset, offsets = 1:length(range(dataset.cell_key)))
    models = map(offsets) do time_offset
        lattice = get_lattice(dataset, time_offset)
        ShroffCelegansModels.build_celegans_model(lattice)
    end
    smodels = map(models) do model
        if ismissing(model)
            missing
        else
            ShroffCelegansModels.Types.StraightenedCelegansModel(model)
        end
    end
    return models, smodels
end

models, smodels = build_models_over_time(datasets["RW10598"][2]);

function average_sliders(models::Vector{<:AbstractCelegansModel})
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    sliders = SliderGrid(f[2,1],
        (label="Weight", range=0:0.01:1),
        (label="Upsampling", range=0:4)
    )

    model = ShroffCelegansModels.average(models, AnalyticWeights([1.0, 0.0]), n_upsample=0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))

    on(throttle(0.1, sliders.sliders[1].value)) do value
        model = ShroffCelegansModels.average(models, AnalyticWeights([1-value, value]); n_upsample=sliders.sliders[2].value[])
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = ShroffCelegansModels.average(models, AnalyticWeights([1-value, value]); n_upsample=upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "weights: ($(round(1-value; digits=2)), $value); number of cross sections: $n_sections"
    end

    display(f)
end

function cross_sections_at_knots(model)
    pts = interpolation_points(model.central_spline)
    sections = map(transverse_splines(model)) do spline
        cs_pts = spline.(pts)
    end
    sections = stack(sections)
    sections = [sections @view(sections[:,1]) fill(Point3(NaN), size(sections,1), 1)]
    sections = PermutedDimsArray(sections, (2,1))
    return vec(sections)
end

function seam_cell_pts(model, n_upsample)
    pts = interpolation_points(model)[1:2^n_upsample:end]
    _transverse_splines = transverse_splines(model)
    [_transverse_splines[1].(pts); _transverse_splines[17].(pts)]
end

#=
function tps_deform_old(x2::AbstractArray,tps::ThinPlateSpline) where T<:Any
	x1,d,c=tps.x1,tps.d,tps.c
	d==[] && throw(ArgumentError("Affine component not available; run tps_solve with compute_affine=true."))

	# deform
	y2 = zeros(eltype(x2),size(x2,1),size(x2,2)+1)
	for i = 1 : size(x2,1)
		z = x2[i:i,:]
		defc = zeros(eltype(x2),1,size(x2,2)+1)
		for j = 1 : size(x1,1)
			n = my_normz - x1[j:j,:])
			defc += tps_basis(n)*c[j:j,:]
		end
		y2[i:i,:] = cat(dims=2, 1.0, z)*d + defc
	end

	y2[:,2:end]
end
=#

function average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    r = range(smts.modelTimeSeries.dataset.cell_key)
    sliders = SliderGrid(f[2,1],
        (label="Time (frames)", range=1.0:0.1:length(r)),
        (label="Upsampling", range=0:4)
    )

    model = smts(1.0, 0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz.(seam_cell_pts(model, 0)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    meshscatter!(ax, _seam_cells; markersize = 10.0, color = :gray, alpha = 0.5)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))

    on(throttle(0.1, sliders.sliders[1].value)) do value
        n_upsample = sliders.sliders[2].value[]
        model = smts(value, n_upsample)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = smts(value, upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, upsampling))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
    end

    display(f)
    f
end

function record_average_sliders(smts::ShroffCelegansModels.StraightenedModelTimeSeries, movieFileName)
    f = average_sliders(smts)
    sg = f.content[2]
    set_close_to!(sg.sliders[2], 2)
    record(f, movieFileName, sg.sliders[1].range[]) do i
        set_close_to!(sg.sliders[1], i)
    end
end

function average_sliders_with_volumes(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    vol_ax = Axis(f[2,1]; xlabel = "Frames", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    r = range(smts.modelTimeSeries.dataset.cell_key)
    sliders = SliderGrid(f[3,1],
        (label="Time (frames)", range=1.0:0.1:length(r)),
        (label="Upsampling", range=0:4)
    )
    slider_range = sliders.sliders[1].range[]

    # precalculate the volumes
    volumes_and_lengths = map(slider_range) do t
        model = smts(t)
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)

    # lengths
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)

    model = smts(1.0, 0)
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz.(seam_cell_pts(model, 0)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    meshscatter!(ax, _seam_cells; markersize = 10.0, color = :gray, alpha = 0.5)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))

    on(throttle(0.1, sliders.sliders[1].value)) do value
        n_upsample = sliders.sliders[2].value[]
        model = smts(value, n_upsample)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
        current_volume[] = Point(value, volumes[floor(Int, (value-1)/step(slider_range) + 1)])
        current_length[] = Point(value, lengths[floor(Int, (value-1)/step(slider_range) + 1)])
    end
    on(throttle(0.1, sliders.sliders[2].value)) do upsampling
        value = sliders.sliders[1].value[]
        model = smts(value, upsampling)
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, upsampling))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
    end

    display(f)
    f
end

function record_average_sliders_with_volumes(smts::ShroffCelegansModels.StraightenedModelTimeSeries, movieFileName)
    f = average_sliders_with_volumes(smts)
    sg = f.content[4]
    set_close_to!(sg.sliders[2], 2)
    record(f, movieFileName, sg.sliders[1].range[]) do i
        set_close_to!(sg.sliders[1], i)
    end
end

normal_flattened_datasets = filter(x->length(range(x.cell_key)) != 259, flattened_datasets)
normal_flattened_datasets = filter(x->x.cell_key.name != "Vab-1_Pos0",normal_flattened_datasets)
smts_datasets = ShroffCelegansModels.StraightenedModelTimeSeries.(normal_flattened_datasets)
lengths = normal_flattened_datasets .|> x->length(range(x.cell_key))
smts_datasets_nt = map(zip(smts_datasets, lengths)) do (ds, _length)
    x -> begin
        nt = x * (_length - 1) + 1.0
        # @info "Normalized time" nt
        ds(nt)
    end
end
models_at_nt(nt) = map(smts_datasets_nt) do ds
    ds(nt)
end
r = LinRange(0.0, 1.0, 201)
avg_models = map(r) do nt
    @info nt
    models = models_at_nt(nt)
    models = filter(!isnothing, models)
    models = identity.(models)
    ShroffCelegansModels.average(models; n_upsample = 2)
end
#=
nothing_count = map(r) do nt
    @info nt
    models = models_at_nt(nt)
    return sum(isnothing.(models))
    models = filter(!isnothing, models)
    models = identity.(models)
    ShroffCelegansModels.average(models; n_upsample = 2)
end
=#
function get_avg_models()
    r = LinRange(0.0, 1.0, 201)
    avg_models = map(r) do nt
        @info nt
        models = models_at_nt(nt)
        models = filter(!isnothing, models)
        models = identity.(models)
        ShroffCelegansModels.average(models; n_upsample = 2)
    end
end

function show_average_models(models)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    vol_ax = Axis(f[2,1]; xlabel = "Time (Normalized)", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    r = LinRange(0.0, 1.0, 201)
    sliders = SliderGrid(f[3,1],
        (label="Time (Normalized)", range=r),
    )
    slider_range = sliders.sliders[1].range[]

    # precalculate the volumes
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)

    # lengths
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)

    n_upsample = 2

    model = models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz))
    _lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))

    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    mesh!(ax, _mesh; colorrange, color = _color, shading)
    meshscatter!(ax, _seam_cells; markersize = 10.0, color = :gray, alpha = 0.5)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 1000))


    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*200 + 1)
        model = models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz)
        _lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        #title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        title[] = "Average over $config_path\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
        current_volume[] = Point(value, volumes[idx])
        current_length[] = Point(value, lengths[idx])
    end

    display(f)
    f
end

function record_show_average_models(models, movieFileName)
    f = show_average_models(models)
    sg = f.content[4]
    record(f, movieFileName, sg.sliders[1].range[]) do i
        set_close_to!(sg.sliders[1], i)
    end
end

const int_cells = ["Mu_int_L", "int1dl", "int1dr", "int1vl", "int1vr", "int2d", "int2v", "int3d", "int3v", "int4d", "int4v", "int5l", "int5r", "int6l", "int6r", "int7l", "int7r", "int8l", "int8r", "int9l", "int9r"]

function get_straightened_annotations(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)::Union{Missing, String}
    timepoint = range(ds.cell_key)[time_offset]
    if timepoint ∈ ds.cell_key.outliers
        return missing
    else
        filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "straightened_annotations", "straightened_annotations.csv")
        if isfile(filepath)
            return filepath
        else
            throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
        end
    end
end

const annotations_cache = Dict{Tuple{String, UnitRange, Bool}, Vector}()

function load_straightened_annotations_over_time(dataset::ShroffCelegansModels.Datasets.NormalizedDataset, offsets::UnitRange = 1:length(range(dataset.cell_key)); use_myuntwist::Bool = false)
    key = (dataset.path, offsets, use_myuntwist)
    if haskey(annotations_cache, key)
        return annotations_cache[key]
    end
    if use_myuntwist
        annotations = map(offsets) do time_offset
            ShroffCelegansModels.untwist_annotations(dataset, time_offset)
        end
        annotations_cache[key] = annotations
        return annotations
    else
        annotations = map(offsets) do time_offset
            path = get_straightened_annotations(dataset, time_offset)
            if ismissing(path)
                return missing
            end
            annotation_df = CSV.read(path, DataFrame)
            pts = Point3f.(eachrow(Matrix(annotation_df)[:, 2:4]))
            pts .-= get_straightened_lattice_xy_center(dataset, time_offset)
            Dict(annotation_df[:,1] .=> pts)
        end
        annotations_cache[key] = annotations
        return annotations
    end
end

function get_model_csv(ds::ShroffCelegansModels.Datasets.NormalizedDataset, data_path, time_offset=1)::Union{Missing, String}
    timepoint = range(ds.cell_key)[time_offset]
    if timepoint ∈ ds.cell_key.outliers
        return missing
    else
        filepath = joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", data_path)
        @info "Filepath" filepath
        if isfile(filepath)
            return filepath
        else
            throw(ArgumentError("$filepath is not a file on disk and is not marked as an outlier."))
        end
    end
end

function get_integrated_annotations(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)
    data_path = joinpath("integrated_annotation","annotations.csv")
    return CSV.read(get_model_csv(ds, data_path, time_offset), DataFrame)
end

function get_integrated_annotations(::Type{Dict}, args...)
    df = get_integrated_annotations(args...)
    return Dict(row[1] => Point3(row[2], row[3], row[4]) for row in eachrow(Matrix(df)))
end

function get_straightened_lattice(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)
    data_path = joinpath("straightened_lattice", "straightened_lattice.csv")
    return CSV.read(get_model_csv(ds, data_path, time_offset), DataFrame)
end

function get_straightened_lattice_xy_center(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset=1)
    csv = get_straightened_lattice(ds, time_offset)
    m = Matrix(csv[1:2:end, 2:3] .+ csv[2:2:end, 2:3])./2
    Point3f(mean(eachrow(m))..., 0)
end

function plot_int_cells(avg_models; flattened_datasets = flattened_datasets)
    int_ds = filter(flattened_datasets) do ds
        "int1dr" in values(ds.cell_key.mapping)
    end

    # models, smodels = build_models_over_time(int_ds[1]);
    
    ds_idx = 1

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(int_ds[ds_idx])
    _length = length(range(int_ds[ds_idx].cell_key))
    smts_nt = x -> begin
        nt = x * (_length - 1) + 1.0
        # @info "Normalized time" nt
        smts(nt, 2)
    end

    smodel = smts_nt(1.0)

    timepoint = 87
    int_annotations_1 = load_straightened_annotations_over_time(int_ds[ds_idx])

    warp_from = ShroffCelegansModels.lattice(smodel) |> vec
    warp_to = ShroffCelegansModels.lattice(avg_models[end]) |> vec
    @info "types" typeof(warp_from) typeof(warp_to)
    tps_solved = tps_solve(warp_from, warp_to, 1)

    pts = collect(values(int_annotations_1[timepoint]))
    @info "pts" typeof(pts)
    transformed_pts = ThinPlateSplines.tps_deform(pts, tps_solved)

    pts = swapyz_scale.(pts)
    transformed_pts = swapyz_scale.(transformed_pts)

    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)

    mesh!(ax, smodel; transform_points = swapyz_scale, alpha = 0.1, transparency = true, color = :green)
    mesh!(ShroffCelegansModels.get_model_contour_mesh(avg_models[timepoint], transform_points = swapyz_scale), alpha = 0.1, transparency = true, color = :magenta)

    scatter!(ax, pts, color = :green)
    scatter!(ax, transformed_pts, color = :magenta)
    
    text!(ax, pts; text = collect(keys(int_annotations_1[timepoint])))
    text!(ax, transformed_pts; text = collect(keys(int_annotations_1[timepoint])))

    arrows!(ax, pts, transformed_pts .- pts)

    ylims!(ax, (0, 1000*voxel_size))
    return f
end

function get_cell_trajectory_dict(dataset::ShroffCelegansModels.Datasets.NormalizedDataset; use_myuntwist = false)
    annotations = load_straightened_annotations_over_time(dataset; use_myuntwist)
    cells = collect(keys(annotations[1]))

    cell_trajectory_dict = map(cells) do cell
        trajectory = map(eachindex(annotations)) do i
            if ismissing(annotations[i]) || !haskey(annotations[i], cell) || any(!isfinite, annotations[i][cell])
                missing
            else
                annotations[i][cell]
            end
        end
        pairs = Iterators.filter(((i,p),)->!ismissing(p),zip(eachindex(trajectory), trajectory))
        pairs = collect(pairs)
        # cell => (first.(pairs), last.(pairs))
        try
            if length(pairs) < 2
                cell => missing
            else
                cell => interpolate((first.(pairs) .- 1)/(length(trajectory)-1), last.(pairs), BSplineKit.BSplineOrder(2))
            end
        catch err
            @error "There was an issue getting cell trajectory" cell dataset.path trajectory
            rethrow()
        end
    end |> Dict

    return cell_trajectory_dict
end

function transform_annotations(from_model, to_model, annotations::AbstractVector{<: Point})
    warp_from = ShroffCelegansModels.lattice(from_model) |> vec
    warp_to = ShroffCelegansModels.lattice(to_model) |> vec
    tps_solved = tps_solve(warp_from, warp_to, 1)

    transformed_pts = ThinPlateSplines.tps_deform(annotations, tps_solved)
    return transformed_pts
end
function transform_annotations(from_model, to_model, annotations::Dict)
    return Dict(keys(annotations) .=> transform_annotations(from_model, to_model, collect(values(annotations))))
end

const annotation_position_cache = Dict{String, Any}()
const my_annotation_position_cache = Dict{String, Any}()

function show_average_models_with_annotations(
    models,
    dataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    #vol_ax = Axis(f[2,1]; xlabel = "Time (Normalized)", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    #length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[2,1],
        (label="Time (Normalized)", range=r),
    )
    slider_range = sliders.sliders[1].range[]

    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
        end
    end

    annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#


    model = models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)

    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = models[idx]
        @info "annotation positions" nt
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end

    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    # _annotation_cells = Observable(annotation_positions(0.0))

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    _annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))


    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        #_lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        #title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        #current_volume[] = Point(value, volumes[idx])
        #current_length[] = Point(value, lengths[idx])
        #smodel = smts_nt(value)
        #=
        _annotation_cells[] = swapyz_scale.(transform_annotations(
            smodel, model, map(values(annotation_dict)) do ann
                ann(value)
            end
        ))
        =#
        _annotation_cells[] = _annotation_positions_over_time[idx]
    end

    display(f)
    f
end

function record_show_average_models_with_annotations(
    models,
    dataset,
    movieFileName;
    use_myuntwist = false
)
    f = show_average_models_with_annotations(models, dataset; use_myuntwist)
    sg = f.content[2]
    record(f, movieFileName, sg.sliders[1].range[]; framerate = 10) do i
        set_close_to!(sg.sliders[1], i)
    end
end

#=
function show_average_models_with_annotations_demo(models, _annotation_positions_over_time)
    f = Figure(size = (1920, 1080))
    title = Observable("Title")
    ax = Axis3(f[1,1], aspect = (1, 10, 1); title)
    #vol_ax = Axis(f[2,1]; xlabel = "Time (Normalized)", ylabel = "Volume (μm^3)", ylabelcolor = :blue)
    #length_ax = Axis(f[2,1]; ylabel = "Length (μm)", yaxisposition = :right, ylabelcolor = :green)
    N_timepoints = 5
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[2,1],
        (label="Time (Normalized)", range=r),
    )
    slider_range = sliders.sliders[1].range[]

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#

    n_upsample = 2

    model = models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(10,0,0)))

    smodel = smts_nt(0.0)

    #=
    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = models[idx]
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end
    =#

    # _annotation_positions_over_time = annotation_positions.(r)

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    # _annotation_cells = Observable(annotation_positions(0.0))


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = collect(keys(annotation_dict)), align = (:right, :bottom))
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))


    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        #_lines[] = swapyz.(cross_sections_at_knots(model))
        _color[] = repeat(color, length(model))
        #title[] = "$(smts.modelTimeSeries.dataset.path)\nt = $value; number of cross sections: $n_sections"
        title[] = "Average over $config_path\nt = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(10,0,0))
        #current_volume[] = Point(value, volumes[idx])
        #current_length[] = Point(value, lengths[idx])
        #smodel = smts_nt(value)
        #=
        _annotation_cells[] = swapyz_scale.(transform_annotations(
            smodel, model, map(values(annotation_dict)) do ann
                ann(value)
            end
        ))
        =#
        _annotation_cells[] = _annotation_positions_over_time[idx]
    end

    display(f)
    f
end
=#
    int_ds = filter(flattened_datasets) do ds
        "int1dr" in values(ds.cell_key.mapping)
    end
    smts2 = ShroffCelegansModels.StraightenedModelTimeSeries(int_ds[2])
    _length2 = length(range(int_ds[2].cell_key))
    smts_nt2 = x -> begin
        nt = x * (_length2 - 1) + 1.0
        # @info "Normalized time" nt
        smts2(nt, 2)
    end

struct NormalizedTimeFunction{F}
    func::F # callable like F
    dataset::ShroffCelegansModels.NormalizedDataset # Parameterize this?
    _length_m_1::Int
    function NormalizedTimeFunction(func::F, dataset) where F
        _length_m_1 = length(range(dataset.cell_key))
        new{F}(func, dataset, _length_m_1)
    end
end
(ntf::NormalizedTimeFunction{F})(nt::Float64) where F = ntf.func(nt * _length_m_1 + 1.0)

"""
    normalize_time(f::Function, dataset::ShroffCelegansModels.NormalizedDataset)

Convert an arbitrary single argument function to one that takes a normalized
time argument ranging from 0.0 to 1.0. The normalized time is scaled to 1.0 to
`length(range(dataset.cell_key))`.
"""
function normalize_time(f, dataset::ShroffCelegansModels.NormalizedDataset)
    return NormalizedTimeFunction(f, dataset)
end
function normalize_time(smts::ShroffCelegansModels.StraightenedModelTimeSeries)
    return normalize_time(smts, smts.modelTimeSeries.dataset)
end

function mipav_df_to_points(df::DataFrame)
    map(df.x_voxels, df.y_voxels, df.z_voxels) do x,y,z
        Point3(x,y,z)
    end
end

function mipav_df_to_point_dict(df::DataFrame)
    map(df.name, df.x_voxels, df.y_voxels, df.z_voxels) do name, x,y,z
        name => Point3(x,y,z)
    end |> Dict
end

second(x) = x[2]

function debug_average_models_with_annotations(
    avg_models,
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset;
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    f = Figure(size = (1920, 1080))
    
    title_prewarp = Observable("Prewarp")
    ax_prewarp = Axis3(f[1,2], aspect = (1, 10, 1); title = title_prewarp)
    
    title = Observable("Title")
    ax = Axis3(f[2,2], aspect = (1, 10, 1); title)

    title_twisted = Observable("Twisted")
    ax_twisted = Axis3(f[1:2,1]; title = title_twisted)

    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[3,1:2],
        (label="Time (Normalized)", range=r),
        (label="Exp. Factor", range=1.0:0.01:1.5),
    )
    slider_range = sliders.sliders[1].range[]

    annotation_text_toggle = Toggle(f)
    central_spline_lines_toggle = Toggle(f)
    f[4, 1:2] = grid!(
        [Label(f, "Annnotation text") annotation_text_toggle Label(f, "Central spline lines") central_spline_lines_toggle]
    )



    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
        end
    end

    mts = smts.modelTimeSeries
    mts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            nt = round(Int, nt)
            title_twisted[] = "Twisted; idx = $nt"
            mts(nt)
        end
    end


    annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#

    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))
    twisted_central_spline = Observable(swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel)))))


    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = avg_models[idx]
        @info "annotation positions" nt
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            ann(nt)
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    twisted_positions = twisted_annotation_positions(0.0)
    twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    # @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    expansion_factor = sliders.sliders[2].value

    twisted_central_pts = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
    twisted_central_pts = Observable(twisted_central_pts)
    # @info twisted_central_pts[]

    twisted_central_line_match = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
    twisted_central_line_match = Observable(twisted_central_line_match)

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    straight_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    twisted_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions))))
    twisted_annotation_text = Observable(twisted_annotation_text)


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    straight_color = Observable(repeat(color, length(smodel)))
    twisted_color = Observable(repeat(color, length(tmodel)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    ann_txt = text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))

    mesh!(ax_prewarp, straight_mesh; colorrange, color = straight_color, transparency = true, alpha = 0.1)
    meshscatter!(ax_prewarp, straight_seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax_prewarp, straight_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax_prewarp, straight_seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    ann_txt = text!(ax_prewarp, straight_annotation_cells; text = straight_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    #lines!(ax_prewarp, _lines, color = :black)
    ylims!(ax_prewarp, (0, 200))

    twisted_seam_cell_text = Observable(String.([tmodel.names[1:2:end]; tmodel.names[2:2:end]]))

    mesh!(ax_twisted, twisted_mesh; colorrange, color = twisted_color, transparency = true, alpha = 0.5)
    lines!(ax_twisted, twisted_central_spline)
    scatter!(ax_twisted, twisted_central_pts)
    cs_lines = lines!(ax_twisted, twisted_central_line_match)
    connect!(cs_lines.visible, central_spline_lines_toggle.active)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        _annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]

        tmodel = mts_nt(value)
        if !ismissing(tmodel)
            twisted_mesh[] = ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale)
            twisted_central_spline[] = swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel))))
            twisted_color[] = repeat(color, length(tmodel))
            twisted_seam_cells[] = swapyz_scale.(seam_cell_pts(tmodel, 0))
            twisted_seam_cell_labels[] = twisted_seam_cells[] .- Ref(Point3f(2,0,0))
            twisted_seam_cell_text[] = [String.(tmodel.names[1:2:end]); String.(tmodel.names[2:2:end])]
            annotation_positions = twisted_annotation_positions(value)
            if !ismissing(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
            end
        end
    end

    on(throttle(0.1, sliders.sliders[2].value)) do expansion_factor_value
        value = sliders.sliders[1].value[]
        tmodel = mts_nt(value)
        annotation_positions = twisted_annotation_positions(value)
        if !ismissing(tmodel)
            if !ismissing(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor_value)))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
            end
        end
    end

    display(f)
    f
end

function debug_average_models_with_annotations(
    avg_models,
    datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset};
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)
    dataset = datasets[1]

    f = Figure(size = (1920, 1080))
    
    title_prewarp = Observable("Prewarp")
    ax_prewarp = Axis3(f[1,2], aspect = (1, 10, 1); title = title_prewarp)
    
    title = Observable("Title")
    ax = Axis3(f[2,2], aspect = (1, 10, 1); title)

    title_twisted = Observable("Twisted")
    ax_twisted = Axis3(f[1:2,1]; title = title_twisted)

    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[3,1:2],
        (label="Time (Normalized)", range=r),
        (label="Exp. Factor", range=1.0:0.01:1.5)
    )
    slider_range = sliders.sliders[1].range[]

    dataset_paths = map(x->x.path, datasets)
    dataset_menu = Menu(f[4,1:2], options = zip(dataset_paths, datasets), default = datasets[2].path)

    dataset = dataset_menu.selection[]

    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
        end
    end

    mts = smts.modelTimeSeries
    mts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            nt = round(Int, nt)
            title_twisted[] = "Twisted; idx = $nt"
            mts(nt)
        end
    end


    annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

    # precalculate the volumes
    #=
    volumes_and_lengths = map(models) do model
        ShroffCelegansModels.volume_by_cross_section(model), ShroffCelegansModels.central_spline(model)(1.0)[3]
    end
    volumes = first.(volumes_and_lengths)
    lines!(vol_ax, slider_range, volumes)
    whole_idx = floor.(Int, 1:1/step(slider_range):length(slider_range))
    scatter!(vol_ax, slider_range[whole_idx], volumes[whole_idx])

    current_volume = Observable(Point(sliders.sliders[1].value[], first(volumes)))
    scatter!(vol_ax, current_volume, marker = '∘', color = :red, markersize = 50)
    =#

    # lengths
    #=
    lengths = last.(volumes_and_lengths)
    lines!(length_ax, slider_range, lengths, color = :green)
    scatter!(length_ax, slider_range[whole_idx], lengths[whole_idx], color = :green)

    current_length = Observable(Point(sliders.sliders[1].value[], first(lengths)))
    scatter!(length_ax, current_length, marker='∘', color = :red, markersize=50)
    =#

    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))

    function annotation_positions(nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = avg_models[idx]
        @info "annotation positions" nt
        swapyz_scale.(transform_annotations(
            _smodel, _model, map(values(annotation_dict)) do ann
                ann(nt)
            end
        ))
    end

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            ann(nt)
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    _annotation_cells = Observable(_annotation_positions_over_time[1])
    straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    twisted_positions = twisted_annotation_positions(0.0)
    twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    _annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))
    straight_annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict)))
    twisted_annotation_text = get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions)))


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    straight_color = Observable(repeat(color, length(smodel)))
    twisted_color = Observable(repeat(color, length(tmodel)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax, _seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))

    mesh!(ax_prewarp, straight_mesh; colorrange, color = straight_color, transparency = true, alpha = 0.1)
    meshscatter!(ax_prewarp, straight_seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax_prewarp, straight_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 1)
    text!(ax_prewarp, straight_seam_cell_labels; text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')], align = (:right, :bottom))
    text!(ax_prewarp, straight_annotation_cells; text = straight_annotation_text, align = (:right, :bottom))
    #lines!(ax_prewarp, _lines, color = :black)
    ylims!(ax_prewarp, (0, 200))

    twisted_seam_cell_text = Observable([tmodel.names[1:2:end]; tmodel.names[2:2:end]])

    mesh!(ax_twisted, twisted_mesh; colorrange, color = twisted_color, transparency = true, alpha = 0.5)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    on(dataset_menu.selection) do ds
        dataset = ds
        smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
        smts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                smts(nt, 2)
            end
        end

        mts = smts.modelTimeSeries
        mts_nt = let _length = length(range(dataset.cell_key))
            x -> begin
                nt = x * (_length - 1) + 1.0
                nt = round(Int, nt)
                title_twisted[] = "Twisted; idx = $nt"
                mts(nt)
            end
        end


        annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

        twisted_annotation_positions = let _length = length(range(dataset.cell_key))
        function twisted_annotation_positions(nt)
            idx = round(Int, nt*(_length - 1) + 1)
            dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
            Dict(keys(dict) .=> swapyz_scale.(values(dict)))
        end
    end


    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end

    straight_annotation_positions_over_time = straight_annotation_positions.(r)

    end

    on(throttle(0.1, sliders.sliders[1].value)) do value
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        _annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]

        tmodel = mts_nt(value)
        if !ismissing(tmodel)
            twisted_mesh[] = ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale)
            twisted_color[] = repeat(color, length(tmodel))
            twisted_seam_cells[] = swapyz_scale.(seam_cell_pts(tmodel, 0))
            twisted_seam_cell_labels[] = twisted_seam_cells[] .- Ref(Point3f(2,0,0))
            twisted_seam_cell_text[] = [tmodel.names[1:2:end]; tmodel.names[2:2:end]]
            annotation_positions = twisted_annotation_positions(value)
            if !ismissing(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
            end
        end
    end

    display(f)
    f
end

function check_dataset(dataset; show_data_frames = false)
    r = range(dataset.cell_key)
    for t in r
        filepath = ShroffCelegansModels.get_lattice_filepath(dataset, t)
        if !isfile(filepath)
            if t ∉ dataset.cell_key.outliers
                @error "$filepath does not exist and is not an outlier"
                println():wait
            end
            continue
        end
        try
            df = CSV.read(filepath, DataFrame)
            _names = df[:,1]
            repeated_names = String[]
            if !allunique(df, 1)
                prev_name = ""
                for name in sort(_names)
                    if name == prev_name
                        push!(repeated_names, name)
                    end
                    prev_name = name
                end
            end
            missing_names = String[]
            for psc in ShroffCelegansModels.psuedo_seam_cells
                if psc ∉ _names
                    push!(missing_names, psc)
                end
                psc_R = replace(psc, 'L' => 'R')
                if psc_R ∉ _names
                    push!(missing_names, psc_R)
                end
            end
            if !isempty(repeated_names) || !isempty(missing_names)
                if !isempty(repeated_names) && !isempty(missing_names)
                    @error "Name error in $filepath" repeated_names missing_names
                elseif !isempty(repeated_names)
                    @error "Repeated names in $filepath" repeated_names
                elseif !isempty(missing_names)
                    @error "Missing names in $filepath" missing_names
                end
                show_data_frames && display(df)
                println()
            end
        catch err
            @error "Error parsing $filepath" err
        end
    end
end

using HDF5

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

using Printf

function save_spline_interpolation(
    parent::Union{HDF5.File, HDF5.Group},
    name::String,
    spline::SplineInterpolation,
    description::String = ""
)
    pts = BSplineKit.SplineInterpolations.interpolation_points(model.central_spline)
    collocation_points = spline.(pts)
    order = BSplineKit.order(spline)
    g = create_group(parent, name)
    g["abscissa"] = pts
    g["ordinate"] = collocation_points
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

#=
function HDF5.hdf5_type_id(::Type{String3})
    HDF5.API.h5t_create(HDF5.API.H5T_STRING,4);
end
function HDF5.datatype(t::Type{String3})
    return HDF5.Datatype(HDF5.hdf5_type_id(t))
end
=#