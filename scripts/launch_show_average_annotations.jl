using Pkg
using Revise
cd(dirname(@__DIR__))
Pkg.activate(dirname(@__DIR__))
using JSON3
using ShroffCelegansModels
using Printf
using HDF5

using InteractiveUtils

@info "Loading demo_averaging.jl..."
@time_imports include("../src/demo_averaging/read_config_json.jl")
@time_imports include("../src/demo_averaging/modelio.jl")
# @time_imports include("../src/demo_averaging.jl")
@info "Loading data..."
@time_imports include("../src/demo_averaging/loading.jl")
@info "Loading show_average_annotations.jl..."
@time_imports include("../src/demo_averaging/seam_cell_pts.jl")
@time_imports include("../src/demo_averaging/load_straightened_annotations_over_time.jl")
@time_imports include("../src/demo_averaging/get_cell_trajectory_dict.jl")
@time_imports include("../src/demo_averaging/show_average_annotations.jl")

const keep_running = Ref(true)

function select_dataset()
    fig = Figure()
    label = Makie.Label(fig[1,1], "Please select a dataset:")
    menu = Menu(fig[2,1], options = collect(keys(datasets)))
    on(menu.selection) do selection
        @info "Launching show_average_annotations(...)" selection
        show_average_annotations(avg_models, datasets[selection]; use_myuntwist=true);
    end
    button = Makie.Button(fig[3,1], label = "Quit")
    on(button.clicks) do b
        @info "Qutting..."
        global keep_running
        keep_running[] = false
        GLMakie.closeall()
        return nothing
    end
    @info "Displaying menu."
    return display(fig)
end

while keep_running[]
    ds = select_dataset()
    if !isnothing(ds)
        wait(ds)
    end
end

println()
@info "Press any key to quit"
readline()

#@info "Launching show_average_annotations(...)"
#show_average_annotations(avg_models, datasets["RW10742"]; use_myuntwist=true);