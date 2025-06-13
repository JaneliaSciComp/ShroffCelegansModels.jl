using Makie
using ShroffCelegansModels
using Printf
using GeometryBasics
using JSON3
using Dates
using Sockets

const ANNOTATION_PERSIST_SERVER_PORT = 3129

function fix_annotation_ap_axis(
    avg_models,
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset;
    use_myuntwist::Bool = false,
    cache::Dict{String} = use_myuntwist ? my_annotation_position_cache : annotation_position_cache,
    ip_address::Sockets.IPAddr = Sockets.getipaddr()
)
    second(x) = x[2]

    f = Figure(size = (1920, 1080))

    title = Label(f[0, 1:2], dataset.path)
    
    title_twisted = Observable("Twisted")
    lscene = LScene(f[1:3,1])
    ax_twisted = lscene.scene

    selected_annotation_name = Observable("")
    ax_distance = Axis(f[1,2],
        xlabel = "Arc distance along central spline",
        ylabel = "Distance to annotation",
        title = selected_annotation_name,
        limits = ((0,200), nothing)
    )
    ax_ratio = Axis(f[2,2],
        xlabel = "Arc distance along central spline",
        ylabel = "Ratio of distance to contour distance",
        title = selected_annotation_name,
        limits = ((0,200), nothing)
    )

    cell_key_range = range(dataset.cell_key)

    ax_z = Axis(f[3,2],
        xlabel = "Time (From Twitch to Hatch)",
        ylabel = "Z",
        title = "Z position of annotation",
        limits = ((first(cell_key_range),last(cell_key_range)), nothing)
    )

    slider_range = LinRange(0.0, 1.0, length(cell_key_range))
    sliders = SliderGrid(f[4, 1:2],
        (label="Time (Normalized)", range=slider_range),
        (label="Timepoint", range=cell_key_range),
        (label="Exp. Factor", range=1.0:0.01:4)
    )
    time_normalized_slider = sliders.sliders[1]
    timepoint_slider = sliders.sliders[2]
    expansion_factor_slider = sliders.sliders[3]

    annotation_text_toggle = Toggle(f)
    central_spline_lines_toggle = Toggle(f)
    contour_mesh_toggle = Toggle(f)
    f[5, 1:2] = grid!(
        [
            Label(f, "Contour mesh");;
            contour_mesh_toggle;;
            Label(f, "Annnotation text");;
            annotation_text_toggle;;
            Label(f, "Central spline lines");;
            central_spline_lines_toggle
        ]
    )
    f[7, :] = buttongrid = GridLayout(tellwidth = false)
    prev_button = buttongrid[1, 1] = Button(f, label = "Previous")
    reset_button = buttongrid[1, 2] = Button(f, label = "Reset")
    update_button = buttongrid[1, 3] = Button(f, label = "Update")
    next_button = buttongrid[1, 4] = Button(f, label = "Next")

    n_upsample = 2

    smts = ShroffCelegansModels.StraightenedModelTimeSeries(dataset)
    smts_nt = let _length = length(range(dataset.cell_key))
        x -> begin
            nt = x * (_length - 1) + 1.0
            smts(nt, 2)
        end
    end

    mts = smts.modelTimeSeries
    mts_nt = let _range = range(dataset.cell_key)
        x -> begin
            nt = x * (length(_range) - 1) + 1.0
            nt = round(Int, nt)
            title_twisted[] = "Twisted; idx = $nt, tp = $(_range[nt])"
            mts(nt)
        end
    end


    annotation_dict = get_cell_trajectory_dict(dataset; use_myuntwist)

    #=
    model = avg_models[1]
    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    smodel = smts_nt(0.0)
    straight_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale))
    straight_seam_cells = Observable(swapyz_scale.(seam_cell_pts(smodel, n_upsample)))
    straight_seam_cell_labels = Observable(straight_seam_cells[] .- Ref(Point3f(2,0,0)))
    =#

    tmodel = mts_nt(0.0)
    twisted_mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(tmodel; transform_points=swapyz_scale))
    twisted_seam_cells = Observable(swapyz_scale.(seam_cell_pts(tmodel, 0)))
    twisted_seam_cell_labels = Observable(twisted_seam_cells[] .- Ref(Point3f(2,0,0)))
    twisted_central_spline = Observable(swapyz_scale.(ShroffCelegansModels.central_spline(tmodel).(LinRange(0,1,length(tmodel)))))


    #=
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
    =#

    function straight_annotation_positions(nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            if ismissing(ann)
                Point3f(NaN, NaN, NaN)
            else
                ann(nt)
            end
        end)
    end

    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end


    #=
    if haskey(cache, dataset.path)
        _annotation_positions_over_time = cache[dataset.path]
    else
        _annotation_positions_over_time = annotation_positions.(r)
        cache[dataset.path] = _annotation_positions_over_time
    end
    =#

    straight_annotation_positions_over_time = straight_annotation_positions.(slider_range)
    # Original positions for reset button
    original_annotation_positions_over_time = deepcopy(straight_annotation_positions_over_time)

    #_annotation_cells = Observable(_annotation_positions_over_time[1])
    straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    twisted_positions = twisted_annotation_positions(0.0)
    twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    # @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    expansion_factor = expansion_factor_slider.value

    twisted_central_pts = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
    twisted_central_pts = Observable(twisted_central_pts)
    # @info twisted_central_pts[]

    twisted_central_line_match = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
    twisted_central_line_match = Observable(twisted_central_line_match)

    function get_annotation_text(dict)
        String.(get.((dataset.cell_key.mapping,), Symbol.(keys(dict)), String.(keys(dict))))
    end

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    #_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    #straight_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    #twisted_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions))))
    _annotation_text = get_annotation_text(annotation_dict)
    straight_annotation_text = get_annotation_text(annotation_dict)
    twisted_annotation_text = get_annotation_text(twisted_positions)
    twisted_annotation_text = Observable(twisted_annotation_text)

    # annotation_menu = Menu(f[6, 1:2], options = twisted_annotation_text)
    menu_options = sort!(collect(values(dataset.cell_key.mapping)))
    annotation_menu = Menu(f[6, 1:2], options = menu_options)
    initial_selected_idx = findfirst(==(first(menu_options)), twisted_annotation_text[])
    if isnothing(initial_selected_idx)
        @error "Could not locate first menu option in twisted_annotation_text" first(menu_options)
        initial_selected_idx = 1
    end

    n_ellipse_pts = length(transverse_splines(tmodel))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    twisted_color = Observable(repeat(color, length(tmodel)))

    twisted_seam_cell_text = Observable(String.([tmodel.names[2:2:end]; tmodel.names[1:2:end]]))

    twisted_mesh_plot = mesh!(
        ax_twisted,
        twisted_mesh;
        colorrange,
        color = twisted_color,
        transparency = true,
        alpha = 0.5,
        inspectable = false
    )
    connect!(twisted_mesh_plot.visible, contour_mesh_toggle.active)
    lines!(ax_twisted, twisted_central_spline)
    scatter!(ax_twisted, twisted_central_pts)
    cs_lines = lines!(ax_twisted, twisted_central_line_match)
    connect!(cs_lines.visible, central_spline_lines_toggle.active)
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    ms_annotation_cells = meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    @info "twisted_seam_cell_labels" twisted_seam_cell_labels[] tmodel.names

    distances = Observable(Float64[])
    selected_distance = Observable(0.0)
    distance_central_pts = Observable(Point3f[])
    central_spline_arc_lengths = Observable(Float64[])

    selected_annotation_idx = Observable(initial_selected_idx)
    max_r = ShroffCelegansModels.max_radius_function(tmodel)
    Npts = length(tmodel)
    z = LinRange(0, 1, Npts)
    max_distance = Observable(Float64[])
    lines!(ax_distance, central_spline_arc_lengths, distances)
    hlines!(ax_distance, selected_distance)
    lines!(ax_distance, central_spline_arc_lengths, max_distance)
    scatter!(ax_distance, central_spline_arc_lengths, distances, color = distances, colormap = Reverse(:viridis))
    scatter!(ax_twisted, distance_central_pts, color = distances, colormap = Reverse(:viridis))
    selected_twisted_annotation_cell = Observable([twisted_annotation_cells[][selected_annotation_idx[]]])
    meshscatter!(ax_twisted, selected_twisted_annotation_cell, color = :red, markersize=1.1)

    ratio = @lift try
        $distances ./ $max_distance
    catch err
        ones(size($max_distance))
    end
    lines!(ax_ratio, central_spline_arc_lengths, ratio)
    hlines!(ax_ratio, 1.0, linestyle = :dash)

    nt_obs = Observable(0.0) 

    on(throttle(0.1, nt_obs)) do value
        #idx = round(Int, value*200 + 1)
        #=
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        # title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        #_annotation_cells[] = _annotation_positions_over_time[idx]

        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        # title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]
        =#

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
                twisted_annotation_text[] = get_annotation_text(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
            end
            #max_r = ShroffCelegansModels.max_radius_function(tmodel)
            #Npts = length(tmodel)
            #z = LinRange(0, 1, Npts)
            #plot_distance(1)
            notify(annotation_menu.selection)
        end

    end

    if !isnothing(time_normalized_slider)
        on(time_normalized_slider.value) do value
            nt_obs[] = value
            cell_key_tp = round(Int, value * (length(cell_key_range)-1)) + first(cell_key_range)
            if expansion_factor_slider.value[] != cell_key_tp
                set_close_to!(timepoint_slider, cell_key_tp)
            end
        end
    end

    on(timepoint_slider.value) do value
        local f = first(cell_key_range)
        local e = last(cell_key_range)
        local nt = (value - f) / (e - f)
        if !isnothing(time_normalized_slider) && !isapprox(time_normalized_slider.value[], nt)
            set_close_to!(time_normalized_slider, nt)
        end
        selected_z_position[] = z_positions[][value-first(cell_key_range)+1]
        @info "Timepoint slider" value
    end

    on(throttle(0.1, expansion_factor_slider.value)) do expansion_factor_value
        # TODO compute value from timepoint slider
        value = time_normalized_slider.value[]
        tmodel = mts_nt(value)
        annotation_positions = twisted_annotation_positions(value)
        if !ismissing(tmodel)
            if !ismissing(annotation_positions)
                twisted_annotation_text[] = get_annotation_text(annotation_positions)
                twisted_annotation_cells[] = collect(values(annotation_positions))
                twisted_central_pts[] = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor_value)))
                twisted_central_line_match[] = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
                #max_r = ShroffCelegansModels.max_radius_function(tmodel)
                #Npts = length(tmodel)
                #z = LinRange(0, 1, Npts)
                #max_distance[] = max_r.(z) .* voxel_size .* expansion_factor_value
                plot_distance(selected_annotation_idx[])
            end
        end
    end



    function plot_distance(idx)
        selected_annotation_idx[] = idx
        #println(idx)
        try
            selected_twisted_annotation_cell[] = [twisted_annotation_cells[][selected_annotation_idx[]]]
        catch err
            @error "Could not select twisted annotation cell" selection_annotation_idx[]
        end
        cs = ShroffCelegansModels.central_spline(tmodel)
        Npts = length(tmodel)
        z = LinRange(0, 1, Npts)
        central_pts = swapyz_scale.(cs.(z))
        central_spline_arc_lengths.val = [0; cumsum(norm.(diff(central_pts)))]
        distance_central_pts[] = central_pts

        max_r = ShroffCelegansModels.max_radius_function(tmodel)
        expansion_factor_value = expansion_factor_slider.value[]
        max_distance[] = max_r.(z) .* voxel_size .* expansion_factor_value

        pt = twisted_annotation_cells[][idx]
        distances[] = norm.(central_pts .- pt)
        #autolimits!(ax_distance)
        #ylims!(ax_distance, nothing)
        limits!(ax_distance, (0, 200), (0, maximum(distances[])))
        limits!(ax_ratio, (0, 200), (0, maximum(ratio[])))
        selected_distance[] = norm(twisted_central_pts[][idx] - pt)

        selected_annotation_name[] = twisted_annotation_text[][idx]
    end
    plot_distance(::Nothing) = plot_distance(initial_selected_idx)

    on(events(f).mousebutton, priority=2) do event
        # print("Mouse clicked")
        if event.button == Mouse.left && event.action == Mouse.press
            p, idx = pick(f)
            if p == ms_annotation_cells
                #println(twisted_annotation_text[][idx])
                plot_distance(idx)
            elseif p == z_lines
                set_close_to!(timepoint_slider, cell_key_range[idx])
            end
        end
    end

    z_positions = let idx=1
        Observable((x->x[idx][2]).(straight_annotation_positions_over_time))
    end
    original_z_positions = let idx=1
        Observable((x->x[idx][2]).(original_annotation_positions_over_time))
    end
    lines!(ax_z, cell_key_range, original_z_positions; color = :gray, linestyle = :dash, label = "Original Z positions")
    z_lines = lines!(ax_z, cell_key_range, z_positions)
    vlines!(ax_z, timepoint_slider.value, color = :red)
    DataInspector(ax_z)

    selected_z_position = Observable(z_positions[][timepoint_slider.value[]-first(cell_key_range)+1])
    vlines!(ax_distance, selected_z_position; color = :red, linestyle = :dash)
    vlines!(ax_ratio, selected_z_position; color = :red, linestyle = :dash)

    #common_annotations_text = collect(keys(annotation_dict))
    common_annotations_text = _annotation_text
    @info common_annotations_text

    on(annotation_menu.selection) do selected
        idx = findfirst(==(selected), twisted_annotation_text[])
        #if selected_annotation_idx[] != idx
            plot_distance(idx)
        #end
        if !isnothing(selected)
            idx_common = findfirst(==(selected), common_annotations_text)
            try
                z_positions[] = (x->x[idx_common][2]).(straight_annotation_positions_over_time)
                original_z_positions[] = (x->x[idx_common][2]).(original_annotation_positions_over_time)
                selected_z_position[] = z_positions[][timepoint_slider.value[] - first(cell_key_range) + 1]
                #autolimits!(ax_z)
                ylims!(ax_z)
            catch err
                @error "Could not get z_positions" err
            end
        end
    end

    on(next_button.clicks) do _
        delta_distances = diff(distances[])
        N = length(delta_distances)
        start = findfirst(>(selected_z_position[]), central_spline_arc_lengths[])
        prev_delta = 0
        for i in 1:N
            idx = mod1(start+i, N)
            cur_delta = delta_distances[idx]
            if prev_delta < 0 && cur_delta >=0
                selected_z_position[] = central_spline_arc_lengths[][idx]
                # @info "Selected z position" selected_z_position[] "at idx" idx prev_delta cur_delta prev_delta < 0 cur_delta >= 0
                break
            end
            # @info "Next button click" i, idx, start, prev_delta, cur_delta
            prev_delta = cur_delta
        end
    end

    on(prev_button.clicks) do _
        delta_distances = diff(distances[])
        N = length(delta_distances)
        start = findfirst(>(selected_z_position[]), central_spline_arc_lengths[])-1
        prev_delta = -1 
        for i in 1:N
            idx = mod1(start-i, N)
            cur_delta = delta_distances[idx]
            if prev_delta >= 0 && cur_delta < 0
                selected_z_position[] = central_spline_arc_lengths[][idx]
                # @debug "Selected z position" selected_z_position[] "at idx" idx prev_delta cur_delta prev_delta >= 0 cur_delta < 0
                break
            end
            # @debug "Prev button click" i, idx, start, prev_delta, cur_delta
            prev_delta = cur_delta
        end
    end

    function update_z_position!(value = nothing)
        if !isnothing(value)
            selected_z_position[] = value
        end
        time_idx = timepoint_slider.value[] - first(cell_key_range) + 1
        z_positions[][time_idx] = selected_z_position[]
        selected = annotation_menu.selection[]
        if isnothing(selected)
            @warn "No annotation selected, cannot update z position"
            return
        end
        idx_common = findfirst(==(selected), common_annotations_text)
        pt = straight_annotation_positions_over_time[time_idx][idx_common]
        pt = Point3{Float64}(pt[1], selected_z_position[], pt[3])
        straight_annotation_positions_over_time[time_idx][idx_common] = pt
        notify(z_positions)
        begin
            # Send the change to the server
            change = AnnotationChange(
                dataset_path = dataset.path,
                annotation_name = selected,
                timepoint = timepoint_slider.value[],
                original_position = original_annotation_positions_over_time[time_idx][idx_common],
                new_z_position = selected_z_position[],
                ip_address = UInt64(ip_address)
            )
            @info "Sending annotation change to server" change
            try
                client = Sockets.connect("0.0.0.0", ANNOTATION_PERSIST_SERVER_PORT)
                @info "Connected to server"
                JSON3.write(client, change)
                println(client, "{}")
                Sockets.close(client)
            catch e
                @error "Failed to send annotation change" exception=e
                update_button.buttoncolor = :red
                reset_button.buttoncolor = :red
                title.text[] = "Failed to send annotation change: $e"
                title.color = :red
            end
        end
        println("Updated z position for annotation $selected at timepoint $(timepoint_slider.value[]) to $(selected_z_position[])")
    end

    on(update_button.clicks) do _
        update_z_position!()
    end

    on(reset_button.clicks) do _
        z = original_z_positions[][timepoint_slider.value[] - first(cell_key_range) + 1]
        update_z_position!(z)
    end

    notify(annotation_menu.selection)
    f
end

using Sockets
using HDF5

@kwdef struct AnnotationChange
    ip_address::UInt64
    dataset_path::String
    annotation_name::String
    timepoint::Int
    original_position::Point3d
    new_position::Point3d
end
function AnnotationChange(;
    ip_address::UInt64,
    dataset_path::String,
    annotation_name::String,
    timepoint::Int,
    original_position::Point,
    new_z_position::Float64
)
    return AnnotationChange(
        ip_address,
        dataset_path,
        annotation_name,
        timepoint,
        original_position,
        Point3d(original_position[1], new_z_position, original_position[3])
    )
end
function AnnotationChange(
    ip_address::UInt64,
    dataset_path::String,
    annotation_name::String,
    timepoint::Int,
    original_position::Point,
    new_z_position::Float64
)
    new_position = Point3d(original_position[1], new_z_position, original_position[3])
    return AnnotationChange(
        ip_address,
        dataset_path,
        annotation_name,
        timepoint,
        original_position,
        new_position
    )
end

function fix_annotation_ap_axis_persist_server(; port = ANNOTATION_PERSIST_SERVER_PORT)
    server = Sockets.listen(port)
    @info "Server listening on port $port"
    server_running = true
    while server_running
        client = Sockets.accept(server)
        @info "Client connected"
        try
            # Handle the client connection
            # fix_annotation_ap_axis_persist(client)
            while true
                str = readavailable(client)
                if isempty(str)
                    @info "No data received, closing connection"
                    break
                end
                data = JSON3.read(str)
                @info "Received request from client" data
                if isempty(data)
                    @info "No data received, closing connection"
                    break
                elseif haskey(data, :shutdown)
                    if data[:shutdown]
                        server_running = false
                        @info "Shutdown request received, closing server"
                    end
                    break
                else
                    c = JSON3.read(str, AnnotationChange)
                    @info "Received annotation change request" c
                    group_name = annotation_change_group_name(
                        c.dataset_path,
                        c.timepoint,
                        c.annotation_name
                    )
                    @info "Group name for HDF5: $group_name"
                    h5open("annotation_changes.h5", "cw", swmr=true) do h5f
                        if !haskey(h5f, group_name)
                            create_group(h5f, group_name)
                        end
                        h5g = h5f[group_name]
                        column_names = (
                            "timepoint",
                            "original_position",
                            "new_position",
                            "timestamp",
                            "ip_address"
                        )
                        column_count = sum(column_names) do c
                            haskey(h5g, c)
                        end
                        if column_count == 0
                            @info "Creating columns in HDF5 group $group_name"
                            timepoint_group = create_dataset(
                                h5g, "timepoint", Int, (1,), max_dims = (-1,), chunk=(16,)
                            )
                            original_position_group = create_dataset(
                                h5g, "original_position", Float64, (3,1), max_dims = (3,-1), chunk=(3,16)
                            )
                            new_position_group = create_dataset(
                                h5g, "new_position", Float64, (3,1), max_dims = (3,-1), chunk=(3,16)
                            )
                            timestamp_group = create_dataset(
                                h5g, "timestamp", Float64, (1,), max_dims = (-1,), chunk=(16,)
                            )
                            ip_address_group = create_dataset(
                                h5g, "ip_address", UInt64, (1,), max_dims = (-1,), chunk=(16,)
                            )
                        elseif column_count == length(column_names)
                            @info "Columns already exist in HDF5 group $group_name"
                            timepoint_group = h5g["timepoint"]
                            original_position_group = h5g["original_position"]
                            new_position_group = h5g["new_position"]
                            timestamp_group = h5g["timestamp"]
                            ip_address_group = h5g["ip_address"]
                            HDF5.set_extent_dims(timepoint_group, (size(timepoint_group, 1) + 1,))
                            HDF5.set_extent_dims(original_position_group, (3, size(original_position_group, 2) + 1))
                            HDF5.set_extent_dims(new_position_group, (3, size(new_position_group, 2) + 1))
                            HDF5.set_extent_dims(timestamp_group, (size(timestamp_group, 1) + 1,))
                            HDF5.set_extent_dims(ip_address_group, (size(ip_address_group, 1) + 1,))
                        else
                            @error "HDF5 group $group_name has inconsistent columns"
                        end
                        @info "Writing annotation change to HDF5 group $group_name"
                        timepoint_group[end] = c.timepoint
                        original_position_group[:,end] = Vector(c.original_position)
                        new_position_group[:,end] = Vector(c.new_position)
                        timestamp_group[end] = Dates.datetime2unix(Dates.now())
                        ip_address_group[end] = UInt64(c.ip_address)
                        @info "Annotation change written to HDF5 group $group_name"
                    end
                end
            end
            #break
        catch e
            bt = Base.catch_backtrace()
            @error "Error handling client: $e" exception=(e,bt)
        finally
            Sockets.close(client)
            @info "Client disconnected"
        end
    end
    close(server)
    @info "Server closed"
    return nothing
end

function annotation_change_group_name(dataset_path, timepoint, annotation_name)
    parts = splitpath(dataset_path)
    parts[1] = replace(parts[1], ":" => "", "\\" => "")
    group_name = join(parts, "/")
    group_name *= "/" * @sprintf("%03d", timepoint)
    group_name *= "/" * annotation_name
    return group_name
end

function shutdown_server()
    s = Sockets.connect("0.0.0.0", ANNOTATION_PERSIST_SERVER_PORT)
    JSON3.write(s, Dict(:shutdown => true))
    close(s)
end