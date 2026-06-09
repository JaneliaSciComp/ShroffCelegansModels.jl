using Pkg
#using Revise
cd(dirname(@__DIR__))
#Pkg.activate(dirname(@__DIR__))
using ShroffCelegansModels
using ShroffCelegansModels.JSON3
using ShroffCelegansModels.Printf
using ShroffCelegansModels.HDF5
using ShroffCelegansModels.ThinPlateSplines # tps_solve
using InteractiveUtils

# The package now owns these symbols (via src/ShroffCelegansModels.jl and
# its transitive includes through parse_worm_dataset_path.jl). The
# include() lines previously used to pull them into Main are gone — we
# import what the still-Main-included files (modelio.jl, loading.jl,
# show_average_annotations.jl, debug_annotation_ap_axis.jl) and this
# launch script's own body actually reference.
using ShroffCelegansModels:
    config_path,
    voxel_size,
    read_config_json,
    seam_cell_pts,
    transform_annotations,
    load_straightened_annotations_over_time,
    get_cell_trajectory_dict,
    get_datasets_info,
    get_group_annotation_positions_over_time,
    interpolation_points,
    second,
    fix_annotation_ap_axis,
    fix_annotation_ap_axis_persist_listen,
    fix_annotation_ap_axis_persist_server,
    update_annotations_cache,
    load_annotation_changes_cache,
    annotations_cache,
    my_annotation_position_cache,
    annotation_position_cache,
    load_avg_models,
    save_annotation_cache,
    load_annotation_cache,
    prime_annotation_caches,
    save_annotations_cache,
    load_annotations_cache

@info "Loading demo_averaging.jl..."
@time_imports include("../src/demo_averaging/modelio.jl")
@info "Loading data..."
@time_imports include("../src/demo_averaging/loading.jl")
@info "Loading show_average_annotations.jl..."
@time_imports include("../src/demo_averaging/show_average_annotations.jl")
@time_imports include("../src/demo_averaging/debug_annotation_ap_axis.jl")

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
    # Map Windows-rooted cache keys (X:\…) to the Linux dataset.path. load_annotation_cache
    # reconstructs the X-tree keys via joinpath on Linux, producing an "X:\/shrofflab/…"
    # form; after the prefix + backslash substitution that leaves a doubled slash
    # ("/nearline/shroff//shrofflab/…"), so collapse it to single so the alias matches
    # the runtime dataset.path (e.g. /nearline/shroff/shrofflab/.../RegB).
    #
    # Use get! so we never overwrite a key already present: the recompute cache's
    # Linux-rooted (nearline) entries load directly under their canonical dataset.path,
    # and those fresh values must win over a legacy X-tree alias that maps to the same
    # path. Re-processing an already-aliased key is a no-op (no X:\ / backslash remain),
    # so iterating while inserting is safe.
    for (k,v) in my_annotation_position_cache
        k2 = replace(k, raw"X:\\" => "$(prefix)")
        k2 = replace(k2, "\\" => "/")
        k2 = replace(k2, "//" => "/")
        get!(my_annotation_position_cache, k2, v)
    end
    for (k,v) in annotations_cache
        a, b, c = k
        a = replace(a, raw"X:\\" => "$(prefix)")
        a = replace(a, "\\" => "/")
        a = replace(a, "//" => "/")
        get!(annotations_cache, (a,b,c), v)
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
else
    alias_cache_unix("/nearline/shroff/")
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
