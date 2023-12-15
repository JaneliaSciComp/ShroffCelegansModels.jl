using BSplineKit.SplineInterpolations: interpolation_points
using ShroffCelegansModels
using ShroffCelegansModels.Types: AbstractCelegansModel
using JSON3
using Missings
using StatsBase
includet("makie.jl")
using GLMakie

const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_full.json"

function read_config_json(config_path::AbstractString = config_path)
    config_json = JSON3.read(config_path)
    cell_keys = Dict{String, Vector{ShroffCelegansModels.CellKey}}()
    map(config_json.data.strains) do strain
        cell_keys[strain.name] = map(strain.folderpaths) do folder_path
            cell_key_json = JSON3.read(joinpath(folder_path, "cell_key.json"))
            ShroffCelegansModels.CellKey(cell_key_json)
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

    model = ShroffCelegansModels.average(models, AnalyticWeights([1.0, 0.0]); n_upsample=0)
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
    [model.transverse_splines[1].(pts); model.transverse_splines[17].(pts)]
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
			n = my_norm(z - x1[j:j,:])
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

    model = smts(1.0; n_upsample = 0)
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
        model = smts(value; n_upsample)
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
        model = smts(value; n_upsample = upsampling)
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
