using ProgressMeter

function show_average_annotations(
    avg_models::Vector{<: ShroffCelegansModels.Types.CelegansModel},
    datasets::Vector{ShroffCelegansModels.Datasets.NormalizedDataset};
    use_myuntwist = false,
    cache = use_myuntwist ? my_annotation_position_cache : annotation_position_cache
)

    f = Figure(size = (1920, 1080))
    
    title = Observable("Annotations Averaged Post-Warp")
    ax = Axis3(f[1:3,1], aspect = (1, 10, 1); title)
    polar_ax = PolarAxis(f[1:3,3]; title)

    title_2d_1 = Observable("X")
    ax_2d_1 = Axis(f[1,2]; title = title_2d_1)
    xlims!(ax_2d_1, (0,1))
    #ylims!(ax_2d_1, (-5,5))
    
    title_2d_2 = Observable("Y")
    ax_2d_2 = Axis(f[2,2]; title = title_2d_2)
    xlims!(ax_2d_2, (0,1))
    #ylims!(ax_2d_2, (-5,5))

    title_2d_3 = Observable("Z")
    ax_2d_3 = Axis(f[3,2]; title = title_2d_3, xlabel = "Time (Normalized)")
    xlims!(ax_2d_3, (0,1))
    ylims!(ax_2d_3, (0, 250))
    # xlabel!("Time (Normalized)")


    N_timepoints = 200
    r = LinRange(0.0, 1.0, N_timepoints + 1)
    sliders = SliderGrid(f[4,1:3],
        (label="Time (Normalized)", range=r),
        (label="Smooth r", range=0:0.01:1),
        (label="Smooth θ", range=0:0.01:1),
        (label="Smooth z", range=0:0.01:1),
    )
    slider_range = sliders.sliders[1].range[]


    annotation_text_toggle = Toggle(f)
    pre_warp_toggle = Toggle(f)
    polar_view_toggle = Toggle(f)
    menu_options = Observable(String[""])
    annotation_menu = Menu(f; options = menu_options)
    f[5, 1:3] = grid!(
            [1,1] => Label(f, "Annnotation text"),
            [1,2] => annotation_text_toggle,
            [1,3] => Label(f, "Pre-warp"),
            [1,4] => pre_warp_toggle,
            [1,5] => Label(f, "Annotation"),
            [1,6] => annotation_menu,
            [1,7] => Label(f, "Polar View"),
            [1,8] => polar_view_toggle,
    )



    n_upsample = 2
    # dataset = first(datasets)

    datasets_info = map(datasets) do dataset
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
        _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
        annotation_dict = Dict(_annotation_text .=> values(annotation_dict))
        return (; dataset, smts_nt, mts_nt, annotation_dict)
    end
    common_annotations = intersect(map(datasets_info) do dataset_info
        collect(keys(dataset_info.annotation_dict))
    end...)

    println(common_annotations)
    #annotation_dict = datasets_info[1].annotation_dict

    #=
    # This is not right, average post-warp
    function average_annotation(name, nt)
        positions = map(datasets_info) do dataset_info
            dataset_info.annotation_dict[name](nt)
        end
        return mean(positions)
    end

    average_annotation_dict = Dict(map(common_annotations) do name
        name => nt->average_annotation(name, nt)
    end)
    annotation_dict = average_annotation_dict
    =#

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

    seam_cell_text = [model.names[1:2:end]; replace.(model.names[1:2:end], 'L' => 'R')]
    menu_options[] = [common_annotations; seam_cell_text]

    _mesh = Observable(ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale))
    #_lines = Observable(swapyz.(cross_sections_at_knots(model)))
    _seam_cells = Observable(swapyz_scale.(seam_cell_pts(model, n_upsample)))
    _seam_cell_labels = Observable(_seam_cells[] .- Ref(Point3f(2,0,0)))

    function annotation_positions(smts_nt, annotation_dict, nt)
        _smodel = smts_nt(nt)
        idx = round(Int, nt*N_timepoints + 1)
        _model = avg_models[idx]
        @info "annotation positions" nt
        positions = missing
        try
            positions = swapyz_scale.(transform_annotations(
                _smodel, _model, map(values(annotation_dict)) do ann
                    if ismissing(ann)
                        return Point3(NaN)
                    else
                        ann(nt)
                    end
                end
            ))
        catch err
            @error "A problem occured at $nt with avg_model[$idx]" exception=(err, Base.catch_backtrace())
            positions = [Point3(NaN) for i in eachindex(annotation_dict)]
        end
        return positions
    end

    function straight_annotation_positions(annotation_dict, nt)
        @info "annotation prewarp positions" nt
        swapyz_scale.(map(values(annotation_dict)) do ann
            if ismissing(ann)
                return Point3(NaN)
            else
                ann(nt)
            end
        end)
    end

    #=
    twisted_annotation_positions = let _length = length(range(dataset.cell_key))
    function twisted_annotation_positions(nt)
        idx = round(Int, nt*(_length - 1) + 1)
        dict = ShroffCelegansModels.twisted_annotations(dataset, idx)
        Dict(keys(dict) .=> swapyz_scale.(values(dict)))
    end
    end
    =#

    
    group_annotation_positions_over_time = map(datasets_info) do dataset_info
        dataset = dataset_info.dataset
        annotation_dict = dataset_info.annotation_dict
        smts_nt = dataset_info.smts_nt
        if haskey(cache, dataset.path)
            _annotation_positions_over_time = cache[dataset.path]
        else
            #_annotation_positions_over_time = annotation_positions.((smts_nt,), (annotation_dict,), r)
            _annotation_positions_over_time = Vector{Vector{Point3{Float64}}}(undef, length(r))
            @showprogress Threads.@threads for i in eachindex(r)
                nt = r[i]
                _annotation_positions_over_time[i] = annotation_positions(smts_nt, annotation_dict, nt)
            end
            cache[dataset.path] = _annotation_positions_over_time
        end
        _annotation_positions_over_time # Vector{Vector{Point3{Float64}}}
        map(_annotation_positions_over_time) do positions
            Dict(keys(annotation_dict) .=> positions)
        end
    end # Vector{Vector{Dict{String, Point3{Float64}}}} # dataset, normalized time, name => position

    seam_cell_positions_over_time = map(eachindex(r)) do j
        _model = avg_models[j]
        seam_cell_pts(_model, n_upsample)
    end

    # averaging
    _annotation_positions_over_time = map(eachindex(r)) do j
        map(common_annotations) do name
            mean(map(eachindex(group_annotation_positions_over_time)) do i
                group_annotation_positions_over_time[i][j][name]
            end)
        end
    end
    println("_annotation_positions_over_time")
    typeof(_annotation_positions_over_time) |> println
    #_annotation_positions_over_time = [_annotation_positions_over_time; seam_cell_positions_over_time]
    _smooth_positions_over_time = Observable(_annotation_positions_over_time)

    group_straight_annotation_positions_over_time = map(datasets_info) do dataset_info
        dataset = dataset_info.dataset
        annotation_dict = dataset_info.annotation_dict
        straight_annotation_positions_over_time = straight_annotation_positions.((annotation_dict,), r)
        map(straight_annotation_positions_over_time) do positions
            Dict(keys(annotation_dict) .=> positions)
        end
    end

    # averaging
    straight_annotation_positions_over_time = map(eachindex(r)) do j
        map(common_annotations) do name
            mean(map(eachindex(group_straight_annotation_positions_over_time)) do i
                group_straight_annotation_positions_over_time[i][j][name]
            end)
        end
    end
    _smooth_positions_over_time = Observable(_annotation_positions_over_time)


    _annotation_cells = Observable(_annotation_positions_over_time[1])
    # straight_annotation_cells = Observable(straight_annotation_positions_over_time[1])
    # twisted_positions = twisted_annotation_positions(0.0)
    # twisted_annotation_cells = Observable(collect(values(twisted_positions)))
    # @info twisted_annotation_cells[]
    # _annotation_cells = Observable(annotation_positions(0.0))

    # expansion_factor = sliders.sliders[2].value
    expansion_factor = 2.0

    # twisted_central_pts = swapyz_scale.(second(ShroffCelegansModels.nearest_central_pt(tmodel, swapyz_unscale.(twisted_annotation_cells[]), expansion_factor[])))
    # twisted_central_pts = Observable(twisted_central_pts)
    # @info twisted_central_pts[]

    # twisted_central_line_match = vec(stack([twisted_annotation_cells[], twisted_central_pts[], fill!(similar(twisted_central_pts[]), Point3f(NaN))]; dims=1))
    # twisted_central_line_match = Observable(twisted_central_line_match)

    #_annotation_text = getindex.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)))
    # _annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    _annotation_text = common_annotations
    # straight_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(annotation_dict)), String.(keys(annotation_dict))))
    # twisted_annotation_text = String.(get.((dataset.cell_key.mapping,), Symbol.(keys(twisted_positions)), String.(keys(twisted_positions))))
    # twisted_annotation_text = Observable(twisted_annotation_text)


    n_ellipse_pts = length(transverse_splines(model))
    colorscheme = :cyclic_wrwbw_40_90_c42_n256
    #shading = MakieCore.automatic
    color=colorschemes[colorscheme][1:256÷n_ellipse_pts:256]
    colorrange = (1,n_ellipse_pts)

    _color = Observable(repeat(color, length(model)))
    # straight_color = Observable(repeat(color, length(smodel)))
    # twisted_color = Observable(repeat(color, length(tmodel)))

    #mesh!(ax, _mesh; colorrange, color = _color, shading, transparency = true)
    scatter_color = use_myuntwist ? :gold : :blue
    mesh!(ax, _mesh; colorrange, color = _color, transparency = true, alpha = 0.1)
    meshscatter!(ax, _seam_cells; markersize = 1.0, color = :gray, alpha = 1)
    meshscatter!(ax, _annotation_cells; markersize = 1.0, color = scatter_color, alpha = 1)
    _selected_annotation = Observable(Point3(NaN))
    meshscatter!(ax, _selected_annotation; markersize = 1.0, color = :red, alpha = 1)
    text!(ax, _seam_cell_labels; text = seam_cell_text, align = (:right, :bottom))
    ann_txt = text!(ax, _annotation_cells; text = _annotation_text, align = (:right, :bottom))
    connect!(ann_txt.visible, annotation_text_toggle.active)
    #lines!(ax, _lines, color = :black)
    ylims!(ax, (0, 200))

    #=

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
    =#

    polar_annotations = lift(_annotation_cells) do positions
        pfc = PolarFromCartesian()
        map(positions) do position
            polar_coord = pfc(Point2(position[1],position[3]))
            Point2(polar_coord.θ, polar_coord.r)
        end
    end
    polar_annotation_selected = lift(_selected_annotation) do position
        pfc = PolarFromCartesian()
        polar_coord = pfc(Point2(position[1], position[3]))
        Point2(polar_coord.θ, polar_coord.r)
    end
    scatter!(polar_ax, polar_annotations; color = scatter_color)
    scatter!(polar_ax, polar_annotation_selected; color = :red)
    text!(polar_ax, polar_annotations; text = common_annotations)

    # selected annotation index
    a = nothing
    vlines!(ax_2d_1, sliders.sliders[1].value, color = :grey)
    vlines!(ax_2d_2, sliders.sliders[1].value, color = :grey)
    vlines!(ax_2d_3, sliders.sliders[1].value, color = :grey)

    on(throttle(0.1, sliders.sliders[1].value)) do value
        set_close_to!(sliders.sliders[2], 0.05)
        set_close_to!(sliders.sliders[3], 0.07)
        set_close_to!(sliders.sliders[4], 0.04)
        idx = round(Int, value*N_timepoints + 1)
        model = avg_models[idx]
        n_sections = length(interpolation_points(model.central_spline))
        _mesh[] = ShroffCelegansModels.get_model_contour_mesh(model; transform_points=swapyz_scale)
        _color[] = repeat(color, length(model))
        # title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $n_sections"
        _seam_cells[] = swapyz_scale.(seam_cell_pts(model, n_upsample))
        _seam_cell_labels[] = _seam_cells[] .- Ref(Point3f(2,0,0))
        if sliders.sliders[2].value[] > 0
            _annotation_cells[] = _smooth_positions_over_time[][idx]
            if isnothing(a)
            elseif a <= length(common_annotations)
                _selected_annotation[] = _smooth_positions_over_time[][idx][a]
            else
            end
        else
            _annotation_cells[] = _annotation_positions_over_time[idx]
            if isnothing(a)
            elseif a <= length(common_annotations)
                _selected_annotation[] = _annotation_positions_over_time[idx][a]
            else
            end
        end



        #=
        smodel = smts_nt(value)
        sn_sections = length(interpolation_points(model.central_spline))
        straight_mesh[] = ShroffCelegansModels.get_model_contour_mesh(smodel; transform_points=swapyz_scale)
        straight_color[] = repeat(color, length(smodel))
        title[] = "Average over $config_path\n$(dataset.path), t = $value; number of cross sections: $sn_sections"
        straight_seam_cells[] = swapyz_scale.(seam_cell_pts(smodel, n_upsample))
        straight_seam_cell_labels[] = straight_seam_cells[] .- Ref(Point3f(2,0,0))
        straight_annotation_cells[] = straight_annotation_positions_over_time[idx]
        =#

        #=
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
        =#
    end

    on(sliders.sliders[3].value) do _
        notify(sliders.sliders[2].value)
    end

    on(sliders.sliders[4].value) do _
        notify(sliders.sliders[2].value)
    end

    on(throttle(0.1, sliders.sliders[2].value)) do smooth_factor
        value = sliders.sliders[1].value[]
        idx = round(Int, value*200 + 1)
        smooth_factor_θ = sliders.sliders[3].value[]
        smooth_factor_z = sliders.sliders[4].value[]
        if smooth_factor > 0 || smooth_factor_θ > 0 || smooth_factor_z > 0
            # For each annotation, get a smoothed radial position
            _smooth_positions_by_annotation = map(eachindex(common_annotations)) do ann_idx
                positions = map(_annotation_positions_over_time) do position
                    position[ann_idx]
                end
                #positions = smooth_radial(positions, 1/smooth_factor)
                positions = smooth_polar_dct1(positions, 1/smooth_factor, 1/smooth_factor_θ, 1/smooth_factor_z)
            end
            # Swap the indexing so that the primary index is time, then annotation
            _smooth_positions_over_time[] = map(eachindex(_annotation_positions_over_time)) do idx
                map(_smooth_positions_by_annotation) do position
                    position[idx]
                end
            end
                _annotation_cells[] = _smooth_positions_over_time[][idx]
            if isnothing(a)
            elseif a <= length(common_annotations)
                _selected_annotation[] = _smooth_positions_over_time[][idx][a]
            end
        else
            _annotation_cells[] = _annotation_positions_over_time[idx]
            if isnothing(a)
            elseif a <= length(common_annotations)
                _selected_annotation[] = _annotation_positions_over_time[idx][a]
            end
        end

        notify(annotation_menu.selection)
        #=
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
        =#
    end

    X = Observable(zeros(size(r)))
    lines!(ax_2d_1, r, X, color = :red)
    Y = Observable(zeros(size(r)))
    lines!(ax_2d_2, r, Y, color = :green)
    Z = Observable(zeros(size(r)))
    lines!(ax_2d_3, r, Z, color = :blue)

    dataset_lines = map(datasets) do datasets
        dsX = Observable(zeros(size(r)))
        lines!(ax_2d_1, r, dsX, color = :gray)
        dsY = Observable(zeros(size(r)))
        lines!(ax_2d_2, r, dsY, color = :gray)
        dsZ = Observable(zeros(size(r)))
        lines!(ax_2d_3, r, dsZ, color = :gray)
        (; dsX, dsY, dsZ)
    end

    # For smoothing
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)

    # pre_warp_selected = false
    on(pre_warp_toggle.active) do active
        notify(annotation_menu.selection)
    end

    on(annotation_menu.selection) do selected
        #set_close_to!(sliders.sliders[2], 0.05)
        #set_close_to!(sliders.sliders[3], 0.07)
        #set_close_to!(sliders.sliders[4], 0.04)
        a = findfirst(==(selected), [common_annotations; seam_cell_text])
        pre_warp_selected = pre_warp_toggle.active[]

        positions_over_time = Point3f[]
        if isnothing(a)
            println("Annotation $a not found")
        elseif a <= length(common_annotations)
            positions_over_time = if pre_warp_selected
                map(straight_annotation_positions_over_time) do positions
                    positions[a]
                end
            else
                map(_annotation_positions_over_time) do positions
                    positions[a]
                end
            end
        else
            a2 = a - length(common_annotations)
            positions_over_time = map(seam_cell_positions_over_time) do positions
                swapyz_scale(positions[a2])
            end
        end
        if !isempty(positions_over_time)
            if sliders.sliders[2].value[] > 0 || sliders.sliders[3].value[] > 0 || sliders.sliders[4].value[] > 0
                # positions_over_time = smooth_radial(positions_over_time, 1/sliders.sliders[2].value[])
                positions_over_time = smooth_polar_dct1(positions_over_time, 1/sliders.sliders[2].value[], 1/sliders.sliders[3].value[], 1/sliders.sliders[4].value[])
            end
            polar_view_selected = polar_view_toggle.active[]
            if !polar_view_selected
                X[] = first.(positions_over_time)
                Y[] = last.(positions_over_time)
            else
                pfc = PolarFromCartesian()
                polar_coords = map(positions_over_time) do position
                    pfc(Point2(first(position), last(position)))
                end
                X[] = (x->x.r).(polar_coords)
                Y[] = (x->x.θ).(polar_coords)
            end
            Z[] = (x->x[2]).(positions_over_time)
            ylims!(ax_2d_1)
            ylims!(ax_2d_2)

            try
                #if a <= length(common_annotations)
                begin
                    for g in eachindex(group_annotation_positions_over_time)
                        (; dsX, dsY, dsZ) = dataset_lines[g]
                        if a > length(common_annotations)
                            dsX[] = [NaN]
                            dsY[] = [NaN]
                            dsZ[] = [NaN]
                            continue
                        end
                        _positions_over_time = map(group_annotation_positions_over_time[g]) do v_d
                            v_d[selected]
                        end
                        dsZ[] = (x->x[2]).(_positions_over_time)
                        if !polar_view_selected
                            dsX[] = first.(_positions_over_time)
                            dsY[] = last.(_positions_over_time)
                        else
                            pfc = PolarFromCartesian()
                            polar_coords = map(_positions_over_time) do position
                                pfc(Point2(first(position), last(position)))
                            end
                            dsX[] = (x->x.r).(polar_coords)
                            dsY[] = (x->x.θ).(polar_coords)
                        end
                    end

                    if a <= length(common_annotations)

                        value = sliders.sliders[1].value[]
                        idx = round(Int, value*N_timepoints + 1)

                        _selected_annotation[] = if sliders.sliders[2].value[] > 0
                            _smooth_positions_over_time[][idx][a]
                        else
                            _annotation_positions_over_time[idx][a]
                        end

                    end
                end
            catch err
                #println(err)
                rethrow(err)
            end

            if polar_view_selected
                title_2d_1[] = "R ($selected), var = $(var(X[]))"
                title_2d_2[] = "θ ($selected), var = $(var(Y[]))"
                ylims!(ax_2d_2, (-π, π))
                title_2d_3[] = "Z ($selected), var = $(var(Z[]))"
            else
                title_2d_1[] = "X ($selected), var = $(var(X[]))"
                title_2d_2[] = "Y ($selected), var = $(var(Y[]))"
                #ylims!(ax_2d_2, (-5, 5))
                ylims!(ax_2d_1)
                ylims!(ax_2d_2)
                title_2d_3[] = "Z ($selected), var = $(var(Z[]))"
            end
        end
    end
    on(polar_view_toggle.active) do polar_view_selected
        notify(annotation_menu.selection)
        #annotation_menu.selection[] = annotation_menu.selection[]
    end

    display(f)
    f
end

function record_show_average_annotations(movieFileName, avg_models, datasets, selection)
    f = show_average_annotations(avg_models, datasets; use_myuntwist=true)
    sg = f.content[6]
    f.content[10].selection[] = selection
    record(f, movieFileName, sg.sliders[1].range[]; framerate = 12) do i
        set_close_to!(sg.sliders[1], i)
    end
end

function smooth_radial(positions_over_time, σ)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()
    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    _filter = G(σ, length(_r))
    _r = abs.(ifft(fft(_r) .* _filter))
    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
end

function smooth_polar(positions_over_time, σ)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()
    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    _filter = G(σ, length(_r))
    _r = abs.(ifft(fft(_r) .* _filter))

    _cs = exp.(_θ .* 1im)
    _cs = ifft(fft(_cs) .* _filter)
    _θ = angle.(_cs)

    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
end

function smooth_polar_dct1(positions_over_time, σ_r, σ_θ, σ_z = 0)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()

    N = length(positions_over_time)

    # DCT Type I Mirroring
    # positions_over_time = [positions_over_time; positions_over_time[end-1:-1:2]]

    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _r = [_r; _r[end-1:-1:2]]
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    if σ_r > 0
        _filter = G(σ_r, length(_r))
        _r = abs.(ifft(fft(_r) .* _filter))
    end
    _r = @view _r[1:N]

    _cs = exp.(_θ .* 1im)
    _cs = [_cs; _cs[end-1:-1:2]]
    if σ_θ > 0
        _filter = G(σ_θ, length(_cs))
        _cs = ifft(fft(_cs) .* _filter)
    end
    _θ = angle.(_cs)

    _z = [_z; _z[end-1:-1:2]]
    if σ_z > 0
        lpz_filter = G(0.5, length(_z))
        lp_z = real.(ifft(fft(_z) .* lpz_filter))
        _filter = G(σ_z, length(_z))
        _z = real.(ifft(fft(_z .- lp_z) .* _filter))
        _z .+= lp_z
    end
    _z = @view _z[1:N]

    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
    # return @view positions_over_time[1:N]
end


"""
dct type II
"""
function smooth_polar_dct2(positions_over_time, σ)
    G(σ,s) = exp.(-fftfreq(s,s).^2 ./2 ./ σ^2)
    pfc = PolarFromCartesian()
    polar_coords = map(positions_over_time) do position
        pfc(Point2(first(position), last(position)))
    end
    _r = (x -> x.r).(polar_coords)
    _θ = (x -> x.θ).(polar_coords)
    _z = (x -> x[2]).(positions_over_time)
    _filter = G(σ, length(_r))
    _r = idct(dct(_r) .* _filter)

    _cs = exp.(_θ .* 1im)
    _cs = idct(dct(_cs) .* _filter)
    _θ = angle.(_cs)

    cfp = CartesianFromPolar()
    positions_over_time = map(_r, _z, _θ) do r, z, θ
        _cartesian = cfp(Polar(r, θ))
        Point3(_cartesian[1], z, _cartesian[2])
    end
    return positions_over_time
end