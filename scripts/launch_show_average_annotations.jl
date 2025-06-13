using Pkg
#using Revise
cd(dirname(@__DIR__))
Pkg.activate(dirname(@__DIR__))
using JSON3
using ShroffCelegansModels
using Printf
using HDF5
using ThinPlateSplines # tps_solve
using InteractiveUtils

@info "Loading demo_averaging.jl..."
@time_imports include("../src/demo_averaging/read_config_json.jl")
@time_imports include("../src/demo_averaging/modelio.jl")
# @time_imports include("../src/demo_averaging.jl")
@info "Loading data..."
@time_imports include("../src/demo_averaging/loading.jl")
@info "Loading show_average_annotations.jl..."
@time_imports include("../src/demo_averaging/seam_cell_pts.jl")
@time_imports include("../src/demo_averaging/transform_annotations.jl")
@time_imports include("../src/demo_averaging/load_straightened_annotations_over_time.jl")
@time_imports include("../src/demo_averaging/get_cell_trajectory_dict.jl")
@time_imports include("../src/demo_averaging/show_average_annotations.jl")
@time_imports include("../src/demo_averaging/debug_annotation_ap_axis.jl")
@time_imports include("../src/demo_averaging/fix_annotation_ap_axis.jl")

function alias_cache(drive_letter)
    if drive_letter == "X"
        return
    end
    for (k,v) in my_annotation_position_cache
        k2 = replace(k, "X:\\" => "$(drive_letter):\\")
        my_annotation_position_cache[k2] = v
    end
    for (k,v) in annotations_cache
        a, b, c = k
        a = replace(a, "X:\\" => "$(drive_letter):\\")
        k2 = (a,b,c)
        annotations_cache[k2] = v
    end
end

function alias_cache_unix(prefix)
    for (k,v) in my_annotation_position_cache
        k2 = replace(k, raw"X:\\" => "$(prefix)")
        k2 = replace(k2, "\\" => "/")
        my_annotation_position_cache[k2] = v
    end
    for (k,v) in annotations_cache
        a, b, c = k
        a = replace(a, raw"X:\\" => "$(prefix)")
        a = replace(a, "\\" => "/")
        k2 = (a,b,c)
        annotations_cache[k2] = v
    end
end

const keep_running = Ref(true)

function select_dataset()
    fig = Figure(; size=(800,600))
    label = Makie.Label(fig[1,1], "Please select a dataset:")
    menu_options = collect(keys(datasets))
    void_option = "[click to select a dataset]"
    pushfirst!(menu_options, void_option)
    menu = Menu(fig[2,1], options = menu_options, width = 400)
    on(menu.selection) do selection
        if selection == void_option
            return
        end
        @info "Launching show_average_annotations(...)" selection
        fig = show_average_annotations(avg_models, datasets[selection]; use_myuntwist=true);
        display(fig)
    end
    button = Makie.Button(fig[3,1], label = "Quit")
    on(button.clicks) do b
        @info "Qutting..."
        global keep_running
        keep_running[] = false
        if @isdefined(GLMakie)
            GLMakie.closeall()
        end
        return nothing
    end
    @info "Displaying menu."
    return display(fig)
end

if gethostname() == "KITTISOPIKULM-2"
    alias_cache("X")
elseif gethostname() == "vm7249"
    # shroff-data.int.janelia.org
    alias_cache_unix("/nearline/shroff")
end

function main()
    while keep_running[]
        ds = select_dataset()
        if !isnothing(ds)
            wait(ds)
        end
    end

    println()
    @info "Press any key to quit"
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#@info "Launching show_average_annotations(...)"
#show_average_annotations(avg_models, datasets["RW10742"]; use_myuntwist=true);
