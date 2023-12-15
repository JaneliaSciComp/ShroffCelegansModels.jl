using Revise, ShroffCelegansModels

includet("makie.jl")

RW10375_datasets = [
    raw"X:\shrofflab\RW10375\Pos0\SPIMB\Reg_Sample\ForTracking\RegB",
    raw"X:\shrofflab\RW10375\Pos1\SPIMB\Reg_Sample\ForTracking\RegB",
    raw"X:\shrofflab\RW10375\Pos2\SPIMB\Reg_Sample\ForTracking\RegB"
]

RW10375_Datasets = ShroffCelegansModels.Datasets.NormalizedDataset.(RW10375_datasets)

RW10375_inputs = map(RW10375_Datasets) do ds
    timepoint = range(ds.cell_key)[38]
    joinpath(ds.path, "Decon_reg_$(timepoint)", "Decon_reg_$(timepoint)_results")
end

RW10375_inputs_lattice = joinpath.(RW10375_inputs, raw"lattice_final\lattice.csv")
models = ShroffCelegansModels.build_celegans_model.(RW10375_inputs_lattice)
smodels = ShroffCelegansModels.Types.StraightenedCelegansModel.(models)

using GLMakie

let f = Figure()
    ax1 = Axis3(f[1,1])
    ax2 = Axis3(f[2,1])
    ax3 = Axis3(f[1,2])
    ax4 = Axis3(f[2,2])
    mesh!(ax1, models[1])
    mesh!(ax2, models[2])
    mesh!(ax3, models[3])
    display(f)
end


avg_smodel = ShroffCelegansModels.average(smodels)

let f = Figure()
    ax1 = Axis3(f[1,1], aspect = (1, 10, 1))
    ax2 = Axis3(f[2,1], aspect = (1, 10, 1))
    ax3 = Axis3(f[1,2], aspect = (1, 10 ,1))
    ax4 = Axis3(f[2,2], aspect = (1, 10, 1))
    mesh!(ax1, smodels[1]; transform_points = swapyz)
    mesh!(ax2, smodels[2]; transform_points = swapyz)
    mesh!(ax3, smodels[3]; transform_points = swapyz)
    mesh!(ax4, avg_smodel; transform_points =  swapyz)
    # lines!(ax4, avg_smodel, swapyz; color = :black)
    display(f)
end