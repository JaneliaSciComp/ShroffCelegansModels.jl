using ShroffCelegansModels
using ShroffCelegansModels.Makie
using ShroffCelegansModels.Printf
using ShroffCelegansModels.GeometryBasics
using ShroffCelegansModels.JSON3
using ShroffCelegansModels.Dates
using ShroffCelegansModels.Sockets
using ShroffCelegansModels.HDF5


const ANNOTATION_PERSIST_SERVER_PORT = 3129

annotation_changes_path() = get(ENV, "ANNOTATION_CHANGES_PATH", joinpath(@__DIR__, "..", "..", "annotation_changes.h5"))

function fix_annotation_ap_axis(
    avg_models,
    dataset::ShroffCelegansModels.Datasets.NormalizedDataset;
    use_myuntwist::Bool = false,
    cache::Dict{String} = use_myuntwist ? my_annotation_position_cache : annotation_position_cache,
    ip_address::Sockets.IPAddr = Sockets.getipaddr(),
    initial_timepoint::Number = 0.0,
    initial_annotation::Union{String, Nothing} = nothing,
    annotation_timepoint_listener::Union{Function, Nothing} = nothing,
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
        #(label="Time (Normalized)", range=slider_range),
        (label="Timepoint", range=cell_key_range),
        (label="Exp. Factor", range=1.0:0.01:4),
        (label="Z (AP axis)", range=0:0.01:200)
    )
    time_normalized_slider = nothing
    slider_index = 1
    #time_normalized_slider = sliders.sliders[slider_index]
    #slider_index += 1
    timepoint_slider = sliders.sliders[slider_index]
    slider_index += 1
    expansion_factor_slider = sliders.sliders[slider_index]
    slider_index += 1
    z_position_slider = sliders.sliders[slider_index]
    slider_index += 1

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
    prev_button = buttongrid[1, 1] = Makie.Makie.Button(f, label = "Previous")
    reset_button = buttongrid[1, 2] = Makie.Button(f, label = "Reset")
    update_button = buttongrid[1, 3] = Makie.Button(f, label = "Update")
    next_button = buttongrid[1, 4] = Makie.Button(f, label = "Next")

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

    # Check if initial annotation is in menu options
    if !isnothing(initial_annotation) && !(initial_annotation in menu_options)
        @error "Initial annotation $initial_annotation not found in menu options"
        initial_annotation = first(menu_options)
    end
    annotation_menu = Menu(f[6, 1:2], options = menu_options, default = initial_annotation)
    @info "initial annotation" initial_annotation menu_options
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
    #connect!(twisted_mesh_plot.visible, contour_mesh_toggle.active)
    twisted_mesh_plot.visible = false
    on(contour_mesh_toggle.active) do v
        twisted_mesh_plot.visible = v
    end
    lines!(ax_twisted, twisted_central_spline)
    scatter!(ax_twisted, twisted_central_pts)
    cs_lines = lines!(ax_twisted, twisted_central_line_match)
    #connect!(cs_lines.visible, central_spline_lines_toggle.active)
    cs_lines.visible = false
    on(central_spline_lines_toggle.active) do v
        cs_lines.visible = v
    end
    meshscatter!(ax_twisted, twisted_seam_cells; markersize = 1.0, color = :gray, alpha = 0.5, transparency = true)
    ms_annotation_cells = meshscatter!(ax_twisted, twisted_annotation_cells; markersize = 1.0, color = use_myuntwist ? :gold : :blue, alpha = 0.5, transparency = true)
    text!(ax_twisted, twisted_seam_cell_labels; text = twisted_seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax_twisted, twisted_annotation_cells; text = twisted_annotation_text, align = (:right, :bottom))
    #connect!(ann_txt.visible, annotation_text_toggle.active)
    ann_txt.visible = false
    on(annotation_text_toggle.active) do v
        ann_txt.visible = v
    end
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
    distances_lines = lines!(ax_distance, central_spline_arc_lengths[], distances[])
    hlines!(ax_distance, selected_distance)
    max_distance_lines = lines!(ax_distance, central_spline_arc_lengths[], max_distance)
    distances_scatter = scatter!(ax_distance, central_spline_arc_lengths[], distances[], color = distances[], colormap = Reverse(:viridis))
    distance_central_pts_scatter = scatter!(ax_twisted, distance_central_pts, color = distances[], colormap = Reverse(:viridis))
    selected_twisted_annotation_cell = Observable([twisted_annotation_cells[][selected_annotation_idx[]]])
    meshscatter!(ax_twisted, selected_twisted_annotation_cell, color = :red, markersize=1.1)


    #=
    ratio = @lift try
        $distances ./ $max_distance
    catch err
        ones(size($max_distance))
    end
    =#
    ratio = Observable(distances[])
    ratio_lines = lines!(ax_ratio, central_spline_arc_lengths[], ratio)
    # hlines!(ax_ratio, 1.0, linestyle = :dash)
    hlines!(ax_ratio, 1.0, linestyle = :solid)

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
        if isnothing(time_normalized_slider)
            nt_obs[] = nt
        end
        selected_z_position[] = z_positions[][value-first(cell_key_range)+1]
        if !isnothing(annotation_timepoint_listener)
            annotation_timepoint_listener(annotation_menu.selection[],value)
        end
        @info "Timepoint slider" value
    end

    on(throttle(0.1, expansion_factor_slider.value)) do expansion_factor_value
        # TODO compute value from timepoint slider
        # value = time_normalized_slider.value[]
        value = nt_obs[]
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
        println("plot_distance")
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
        # central_spline_arc_lengths.val = [0; cumsum(norm.(diff(central_pts)))]
        _central_spline_arc_lengths = [0; cumsum(norm.(diff(central_pts)))]
        central_spline_arc_lengths[] = _central_spline_arc_lengths
        # distance_central_pts[] = central_pts

        max_r = ShroffCelegansModels.max_radius_function(tmodel)
        expansion_factor_value = expansion_factor_slider.value[]
        _max_distance = max_r.(z) .* voxel_size .* expansion_factor_value
        # max_distance[] = max_r.(z) .* voxel_size .* expansion_factor_value
        Makie.update!(max_distance_lines, arg1 = _central_spline_arc_lengths, arg2 = _max_distance)

        pt = twisted_annotation_cells[][idx]
        _distances = norm.(central_pts .- pt)
        distances[] = _distances
        #distances[] = norm.(central_pts .- pt)
        Makie.update!(distances_lines, arg1 = _central_spline_arc_lengths, arg2 = _distances)
        Makie.update!(distances_scatter, arg1 = _central_spline_arc_lengths, arg2 = _distances, color = _distances)
        Makie.update!(distance_central_pts_scatter, arg1 = central_pts, color = _distances)

        _ratio = _distances ./ _max_distance
        Makie.update!(ratio_lines, arg1 = _central_spline_arc_lengths, arg2 = _ratio)

        #autolimits!(ax_distance)
        #ylims!(ax_distance, nothing)
        limits!(ax_distance, (0, 200), (0, maximum(_distances)))
        limits!(ax_ratio, (0, 200), (0, maximum(_ratio)))
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
    lines!(ax_z, cell_key_range, original_z_positions; color = :gray, linestyle = :solid, label = "Original Z positions")

    z_lines = lines!(ax_z, cell_key_range, z_positions)
    vlines!(ax_z, timepoint_slider.value, color = :red)
    DataInspector(ax_z)

    selected_z_position = Observable(z_positions[][timepoint_slider.value[]-first(cell_key_range)+1])
    vlines!(ax_distance, selected_z_position; color = :red, linestyle = :solid)
    vlines!(ax_ratio, selected_z_position; color = :red, linestyle = :solid)

    #common_annotations_text = collect(keys(annotation_dict))
    common_annotations_text = _annotation_text
    @info common_annotations_text

    h5open(annotation_changes_path(), "r", swmr=true) do h5f
        for time_idx in eachindex(cell_key_range)
            timepoint = cell_key_range[time_idx]
            for annotation_name in values(dataset.cell_key.mapping)
                annotation_idx = findfirst(==(annotation_name), common_annotations_text)
                group_name = annotation_change_group_name(dataset.path, timepoint, annotation_name)
                if haskey(h5f, group_name)
                    @info "Loading annotation changes from annotations_changes.h5" group_name annotation_idx annotation_name
                    try
                        new_position = h5f[group_name]["new_position"][:,end]
                        straight_annotation_positions_over_time[time_idx][annotation_idx] = Point3f(
                            new_position[1],
                            new_position[2],
                            new_position[3]
                        )
                        @info "Loaded" new_position
                    catch err
                        println(err)
                    end
                end
            end
        end
    end


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
            if !isnothing(annotation_timepoint_listener)
                annotation_timepoint_listener(selected, timepoint_slider.value[])
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

    update_z_slider::Bool = true
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
        if update_z_slider
            set_close_to!(z_position_slider, selected_z_position[])
        end
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

    on(z_position_slider.value) do value
        update_z_slider = false
        update_z_position!(value)
        update_z_slider = true
    end

    on(update_button.clicks) do _
        update_z_position!()
    end

    on(reset_button.clicks) do _
        z = original_z_positions[][timepoint_slider.value[] - first(cell_key_range) + 1]
        update_z_position!(z)
    end

    notify(annotation_menu.selection)

    if !isnothing(initial_timepoint) && initial_timepoint != 0.0
        set_close_to!(timepoint_slider, initial_timepoint)
    end

    f
end

struct AnnotationChange
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
    new_position::Union{Point, Nothing} = nothing,
    new_z_position::Union{Float64, Nothing} = nothing,
)
    isnothing(new_position) && isnothing(new_z_position) && throw(ArgumentError("One of new_position or new_z_position must be provided"))
    if isnothing(new_position)
        new_position = Point3d(original_position[1], new_z_position, original_position[3])
    end
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
                    h5open(annotation_changes_path(), "cw", swmr=true) do h5f
                        if !haskey(h5f, group_name)
                            create_group(h5f, group_name)
                        end
                        h5g = h5f[group_name]
                        column_names = (
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
                            original_position_group = create_dataset(
                                h5g, "original_position", Float64, ((3,1), (3,-1)), chunk=(3,16)
                            )
                            new_position_group = create_dataset(
                                h5g, "new_position", Float64, ((3,1), (3,-1)), chunk=(3,16)
                            )
                            timestamp_group = create_dataset(
                                h5g, "timestamp", Float64, ((1,), (-1,)), chunk=(16,)
                            )
                            ip_address_group = create_dataset(
                                h5g, "ip_address", UInt64, ((1,), (-1,)), chunk=(16,)
                            )
                        elseif column_count == length(column_names)
                            @info "Columns already exist in HDF5 group $group_name"
                            original_position_group = h5g["original_position"]
                            new_position_group = h5g["new_position"]
                            timestamp_group = h5g["timestamp"]
                            ip_address_group = h5g["ip_address"]
                            HDF5.set_extent_dims(original_position_group, (3, size(original_position_group, 2) + 1))
                            HDF5.set_extent_dims(new_position_group, (3, size(new_position_group, 2) + 1))
                            HDF5.set_extent_dims(timestamp_group, (size(timestamp_group, 1) + 1,))
                            HDF5.set_extent_dims(ip_address_group, (size(ip_address_group, 1) + 1,))
                        else
                            @error "HDF5 group $group_name has inconsistent columns"
                        end
                        @info "Writing annotation change to HDF5 group $group_name"
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

function load_annotation_changes_cache(filepath = annotation_changes_path())
    # changes = Dict{String,Pair{Vector{Point3{Float64}},Vector{Point3{Float64}}}}()
    changes = Dict{String,Pair{Point3{Float64},Point3{Float64}}}()

    h5open(filepath, "r") do h5f
        function _descend(p::Union{HDF5.File,HDF5.Group})
            for k in keys(p)
                _descend(p[k])
            end
        end
        function _descend(d::HDF5.Dataset)
            _name = HDF5.name(d)
            # @info "Processing dataset" _name d[]
            _paths = splitpath(_name)
            popfirst!(_paths)
            idx = tryparse(Int, _paths[end-2])
            if idx === nothing
                idx = parse(Int, _paths[end-3])
            end
            # @info "Paths after popfirst" idx _paths[end-3:end]

            #idx = pop!(_paths)
            #idx = parse(Int, idx)

            # TODO: check if _paths[1] is nearline or a Windows UNC
            #_paths[1] = _paths[1] * ":\\"
            # _path = join(_paths, '/')
            _path = join([_paths[1] * ":", _paths[3:end]...], '\\')

            if _paths[end] != "new_position"
                # @info "Skipping dataset" _path "not a new_position dataset"
                return
            end

            _old_position_path = join([_paths[1:end-1]..., "original_position"], '/')
            old_positions = h5f[_old_position_path][]
            old_pts = Point3{Float64}.(eachcol(old_positions))

            data = d[]
            pts = Point3{Float64}.(eachcol(data))
            @debug "Got points" pts
            #println(_path)
            #println(HDF5.name(d))
            #=
            cache = get!(my_annotation_position_cache, _path) do
                P = parent(d)
                N = count(keys(P)) do k
                    isa(P[k], HDF5.Dataset)
                end
                Vector{Vector{Point3{Float64}}}(undef, N)
            end
            cache[idx] = pts
            =#
            old_pt = first(old_pts)
            new_pt = last(pts)
            # Ignore changes if the old and new points are the same
            if norm(new_pt - old_pt) < 1
                @info "Ignoring change for $_path, old and new points are the same: $(old_pt) => $(new_pt)"
                return
            end
            _path = replace(_path, "nearline" => "X")
            changes[_path] = old_pt => new_pt
        end
        _descend(h5f)
    end
    return changes
end

"""
    update_annotations_cache(annotations_cache, annotations_changes)

Update the `annotations_cache` with the changes in `annotations_changes`.

Parameters
----------
- `annotation_cache`: The current cache of annotation positions.
    Usually ShroffCelegansModels.annotations_cache.
- `annotations_changes`: The changes to be applied to the annotation cache.
    Usually from `load_annotation_changes_cache`.

Returns
----------
- The updated `annotations_cache`.
"""
function update_annotations_cache(
    annotations_cache::Dict{Tuple{String, UnitRange, Bool}, AnnotationsCacheValue},
    annotations_changes::Dict{String, Pair{Point3{Float64},Point3{Float64}}};
    dry_run::Bool = false
)
    annotations_cache_keys = keys(annotations_cache)
    key_map_dict = Dict(first.(annotations_cache_keys) .=> annotations_cache_keys)
    for (change_key, (old_pt, new_pt)) in annotations_changes
        change_key_parts = split(change_key, "\\")
        timepoint = tryparse(Int, change_key_parts[end-3])
        if isnothing(timepoint)
            # Case when annotation name does not contain a slash
            timepoint = parse(Int, change_key_parts[end-2])
            dataset_path = join(change_key_parts[begin:end-3], "\\")
            annotation_name = change_key_parts[end-1]
        else
            # Case when annotation name contains a slash
            # timepoint is change_key_parts[end-3]
            dataset_path = join(change_key_parts[begin:end-4], "\\")
            annotation_name = change_key_parts[end-2] * "/" * change_key_parts[end-1]
        end
        # Could error if dataset_path is not in key_map_dict
        if !haskey(key_map_dict, dataset_path)
            @warn "Dataset path $dataset_path not found in annotations cache keys, skipping change for annotation $annotation_name at timepoint $timepoint"
            continue
        end
        annotations_cache_key = key_map_dict[dataset_path]
        local_dataset_path = dataset_path
        if Sys.isunix()
            local_dataset_path = replace(local_dataset_path, raw"X:\\" => "/nearline/shroff/")
            local_dataset_path = replace(local_dataset_path, "\\" => "/")
        end
        dataset = ShroffCelegansModels.Dataset(local_dataset_path)
        annotation_symbol = findfirst(==(annotation_name), dataset.cell_key.mapping)
        if isnothing(annotation_symbol)
            @warn "Annotation $annotation_name not found in dataset $local_dataset_path, skipping change at timepoint $timepoint"
            continue
        end
        # timepoint above is the actual timepoint, but we need to adjust it to the dataset's range
        timepoint = timepoint - dataset.cell_key.start + 1
        cached_annotations = annotations_cache[annotations_cache_key].annotations
        if ismissing(cached_annotations[timepoint])
            @warn("Annotation $annotation_name at timepoint $timepoint is missing in cache for $annotations_cache_key")
            if !dry_run
                cached_annotations[timepoint] = Dict{String, Point3{Float64}}()
                cached_annotations[timepoint][string(annotation_symbol)] = swapyz_unscale(new_pt)
            end
        end
        # Check old_pt matches the cached point
        try
            if !isapprox(cached_annotations[timepoint][string(annotation_symbol)], swapyz_unscale(old_pt))
                _norm = norm(cached_annotations[timepoint][string(annotation_symbol)] - swapyz_unscale(old_pt))
                @warn """Annotation point mismatch for $annotations_cache_key at timepoint $timepoint: $(cached_annotations[timepoint][string(annotation_symbol)]) != $(swapyz_unscale(old_pt)), norm difference: $_norm.
                This may indicate the cache is out of sync with the changes, or the change does not apply cleanly to the current cache state. Consider reviewing this change and the current cache state to ensure consistency."""
            end
        catch e
            @warn "Error checking annotation point for $annotations_cache_key at timepoint $timepoint with annotation $annotation_symbol: $e"
        end
        if !dry_run
            cached_annotations[timepoint][string(annotation_symbol)] = swapyz_unscale(new_pt)
        end
    end
    return annotations_cache
end
