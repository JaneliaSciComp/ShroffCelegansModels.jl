using Revise, ShroffCelegansModels, GLMakie

includet("makie.jl")

RW10375_datasets = [
    raw"X:\shrofflab\RW10375\Pos0\SPIMB\Reg_Sample\ForTracking\RegB",
    raw"X:\shrofflab\RW10375\Pos1\SPIMB\Reg_Sample\ForTracking\RegB",
    raw"X:\shrofflab\RW10375\Pos2\SPIMB\Reg_Sample\ForTracking\RegB"
]

RW10375_Datasets = ShroffCelegansModels.Datasets.NormalizedDataset.(RW10375_datasets)

function get_lattice(ds::ShroffCelegansModels.Datasets.NormalizedDataset, time_offset)
    timepoint = range(ds.cell_key)[time_offset]
    joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results", "lattice_final", "lattice.csv")
end

function build_straightend_models(datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset}, time_offset)
    models = map(datasets) do ds
        ShroffCelegansModels.build_celegans_model(get_lattice(ds, time_offset))
    end
    smodels = map(models) do model
        ShroffCelegansModels.Types.StraightenedCelegansModel(model)
    end
    avg_smodel = ShroffCelegansModels.average(smodels)
    return smodels, avg_smodel
end

function show_sliders(datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset}; save_movie = false)
    timepoints = map(datasets) do ds
        range(ds.cell_key)
    end
    min_length = minimum(length.(timepoints))
    f = Figure(size = (1920, 1080))
    ax = Vector{Axis3}(undef, 4)
    titles = Observable.(["Pos 0", "Pos 1", "Pos 2", "Average"])
    ax[1] = Axis3(f[1,1], aspect = (1, 10, 1), title = titles[1])
    ax[2] = Axis3(f[1,2], aspect = (1, 10, 1), title = titles[2])
    ax[3] = Axis3(f[2,1], aspect = (1, 10 ,1), title = titles[3])
    ax[4] = Axis3(f[2,2], aspect = (1, 10, 1), title = titles[4])
    sliders = SliderGrid(f[3,1:2],
        (label="Frame", range=1:min_length)
    )
    smodels, avg_smodel = build_straightend_models(datasets, 1)
    all_models = [smodels..., avg_smodel]
    meshes = ShroffCelegansModels.get_model_contour_mesh.(all_models; transform_points = swapyz)
    obs = Observable.(meshes)

    # cs = Dict(ShroffCelegansModels.cross_sections.(all_models).sections)

    # Display settings
    n_ellipse_pts = length(transverse_splines(avg_smodel))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)
    alpha = 0.9

    rendered_meshes = Vector{MakieCore.Mesh}(undef, 4)
    for i in 1:4
        rendered_meshes[i] = mesh!(ax[i], obs[i]; shading, colorrange, alpha, color=repeat(color, length(all_models[i])))
        ylims!(ax[i], 0, 1000)
    end
    display(f)
    function render_frame(i)
        smodels, avg_smodel = build_straightend_models(datasets, i)
        all_models = [smodels..., avg_smodel]
        meshes = ShroffCelegansModels.get_model_contour_mesh.(all_models; transform_points = swapyz)
        for j in 1:4
            obs[j][] = meshes[j]
            rendered_meshes[j].color[] = repeat(color, length(all_models[j]))
        end
        for j in 1:3
            titles[j][] = "Pos $(j-1): time = $(timepoints[j][i])"
        end
        titles[4][] = "Average: time = $i"
        sleep(0.01)
    end
    on(throttle(0.1, sliders.sliders[1].value)) do value
        render_frame(value)
    end
    #=
    if save_movie
        record(render_frame, f, "RW10375_movie.mp4", 1:min_length, framerate = 5)
    else
        foreach(render_frame, 1:min_length)
    end
    =#
end

show_sliders(RW10375_Datasets)