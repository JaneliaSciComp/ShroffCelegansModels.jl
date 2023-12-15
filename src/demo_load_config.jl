using JSON3
using ShroffCelegansModels
using Missings
using StatsBase
includet("makie.jl")
using GLMakie

const config_path = raw"D:\shroff\python_model_building\C-Elegans-Model-Generation\config_full.json"

function read_config_json(config_path::AbstractString = config_path; exclude = nothing)
    config_json = JSON3.read(config_path)
    cell_keys = Dict{String, Vector{ShroffCelegansModels.CellKey}}()
    map(config_json.data.strains) do strain
        cell_keys[strain.name] = map(filter(!=(exclude), strain.folderpaths)) do folder_path
            cell_key_json = JSON3.read(joinpath(folder_path, "cell_key.json"))
            ShroffCelegansModels.CellKey(cell_key_json)
        end
    end

    datasets = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
    map(config_json.data.strains) do strain
        datasets[strain.name] = map(filter(!=(exclude),strain.folderpaths)) do folder_path
            ShroffCelegansModels.NormalizedDataset(joinpath(folder_path, "RegB"))
        end
    end
    return config_json, cell_keys, datasets
end

exclude = raw"X:\shrofflab\Non-model C elegans folders\Untwisting_Paper\DCR4221\063014_lattices_Javier - UTP 1\Decon_reg"
config_json, cell_keys, datasets = read_config_json(; exclude)

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

function interpolate_models_over_time(models::Vector{<: Union{Missing,ShroffCelegansModels.Types.StraightenedCelegansModel}})
    interpolated_models = Vector{ShroffCelegansModels.Types.AbstractCelegansModel}(undef, length(models))
    for i in eachindex(models)
        if ismissing(models[i])
            start_model = models[i-1]
            start_index = i - 1
            end_model = nothing
            end_index = nothing
            for j in i+1:lastindex(models)
                if !ismissing(models[j])
                    end_model = models[j]
                    end_index = j
                    break
                end
            end
            gap_length = end_index - start_index
            for j in start_index+1:end_index-1
                weights = aweights([j-start_index end_index-j])
                interpolated_models[j] = ShroffCelegansModels.average(disallowmissing(smodels[[start_index,end_index]]), weights)
            end
            # TODO: Resume here, interpolate over missing models
            @info "Start and end" start_index end_index
        else
            interpolated_models[i] = models[i]
        end
    end
    return interpolated_models
end

function average_straightened_models(smodels)
    avg_smodel = ShroffCelegansModels.average(smodels)
    return smodels, avg_smodel
end

function get_averaged_models(flattened_datasets = flattened_datasets)
    min_length = map(flattened_datasets) do ds
        length(range(ds.cell_key))
    end |> minimum
    avg_smodels = map(1:min_length) do offset
        @info "Processing" offset min_length
        try
            smodels, avg_smodel = build_straightend_models(flattened_datasets, offset)
            return avg_smodel
        catch err
            @error err
        end
        return nothing
    end
    return avg_smodels
end

avg_smodels = get_averaged_models(flattened_datasets)

filtered_avg_smodels = filter(!isnothing, avg_smodels);

function plot_models(models; plot_lines = false, plot_cross_sections = false)
    f = Figure(size = (1920, 1080))
    ax = Axis3(f[1,1], aspect = (1,10,1))
    sliders = SliderGrid(f[2,1],
        (label="Frame", range=1:length(models))
    )
    mesh_obs = Observable(ShroffCelegansModels.get_model_contour_mesh(first(models); transform_points = swapyz))

    n_ellipse_pts = length(transverse_splines(first(models)))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)
    alpha = 0.9


    ylims!(ax, 0, 1000)
    rendered_mesh = mesh!(ax, mesh_obs; shading, colorrange, alpha, color=repeat(color, length(first(models))))

    model_obs = Observable{ShroffCelegansModels.Types.AbstractCelegansModel}(first(models))
    if plot_lines
        rendered_lines = lines!(ax, model_obs, swapyz; color = :black)
    end
    if plot_cross_sections
        ShroffCelegansModels.cross_sections(first(models)).sections
    end
    on(throttle(0.1, sliders.sliders[1].value)) do value
        mesh_obs[] = ShroffCelegansModels.get_model_contour_mesh(models[value]; transform_points = swapyz)
        rendered_mesh.color[] = repeat(color, length(models[value]))
        model_obs[] = models[value]
    end
    display(f)
    sliders
end

plot_models(filtered_avg_smodels; plot_lines = true)

models, smodels = build_models_over_time(datasets["RW10598"][2]);
smodels_interpolated = interpolate_models_over_time(smodels);


let f = Figure()
    ax = Axis3(f[1,1], aspect = (1, 10, 1))
    mesh!(ax, smodels_interpolated[13]; transform_points = swapyz, color = :red, alpha = 0.25)
    mesh!(ax, smodels_interpolated[14]; transform_points = swapyz, color = :blue, alpha = 0.5)
    mesh!(ax, smodels_interpolated[15]; transform_points = swapyz, color = :green, alpha = 0.25)
    display(f)
end

let f = Figure()
    ax = Vector{Axis3}(undef, 4)
    ax[1] = Axis3(f[1,1], aspect = (1,10,1))
    mesh!(ax[1], smodels_interpolated[13]; transform_points = swapyz, color = :red, alpha = 0.5)
    text!(ax[1], smodels_interpolated[13], swapyz)
    ylims!(ax[1], 0, 1000)

    ax[2] = Axis3(f[1,2], aspect = (1,10,1))
    mesh!(ax[2], smodels_interpolated[15]; transform_points = swapyz, color = :green)
    text!(ax[2], smodels_interpolated[15], swapyz)
    ylims!(ax[2], 0, 1000)

    ax[3] = Axis3(f[2,1], aspect = (1,10,1))
    mesh!(ax[3], smodels_interpolated[14]; transform_points = swapyz, color = :blue)
    text!(ax[3], smodels_interpolated[14], swapyz)
    ylims!(ax[3], 0, 1000)

    ax[4] = Axis3(f[2,2], aspect = (1,10,1))
    mesh!(ax[4], smodels_interpolated[13]; transform_points = swapyz, color = :red, alpha = 0.25)
    mesh!(ax[4], smodels_interpolated[14]; transform_points = swapyz, color = :blue, alpha = 0.5)
    mesh!(ax[4], smodels_interpolated[15]; transform_points = swapyz, color = :green, alpha = 0.25)
    text!(ax[4], smodels_interpolated[14], swapyz)
    ylims!(ax[4], 0, 1000)
    display(f)
end

let f = Figure()
    smodels = (flattened_datasets .|> x->last(build_models_over_time(x, 1:1))) .|> last;
    avg_smodel1 = ShroffCelegansModels.average(smodels; n_upsample=2)
    smodel1 = ShroffCelegansModels.average([smodels[1]]; n_upsample=2)
    start_pts = stack(filter(!isnan, cross_sections_at_knots(smodel1)); dims=1)
    end_pts = stack(filter(!isnan, cross_sections_at_knots(avg_smodel1)); dims=1)
    tps =  tps_solve(start_pts, end_pts, 1.0; compute_affine=true)

    pts = pts = stack(CartesianIndices((-75:25:75, -75:25:75, 1:100:1000)) .|> x->Point3(x.I); dims=1)
    deformed_pts = tps_deform(pts, tps)

    ax = Vector{Axis3}(undef, 4)
    ax[1] = Axis3(f[1,1], aspect = (1,10,1))
    mesh!(ax[1], smodel1; transform_points = swapyz)
    scatter!(ax[1], pts[:, [1,3,2]])
    ylims!(ax[1], 0, 1000)
    display(f)


    ax[2] = Axis3(f[2,1], aspect = (1,10,1))
    mesh!(ax[2], avg_smodel1; transform_points = swapyz)
    scatter!(ax[2], deformed_pts[:, [1,3,2]]; color = :gold)
    ylims!(ax[2], 0, 1000)
    display(f)

    ax[3] = Axis3(f[3,1], aspect = (1,10,1))
    #mesh!(ax[3], avg_smodel1; transform_points = swapyz)
    scatter!(ax[3], pts[:, [1,3,2]])
    scatter!(ax[3], deformed_pts[:, [1,3,2]]; color = :gold)
    arrows!(ax[3], swapyz.(Point3.(eachrow(pts))), swapyz.(Point3.(eachrow(deformed_pts - pts))); arrowsize)
    ylims!(ax[3], 0, 1000)
    display(f)
end