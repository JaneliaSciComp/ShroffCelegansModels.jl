using Dates

function parse_worm_dataset_path(s::String)
    date_pattern = r"\d{6}"
    position_pattern = r"Pos\d+"
    date_match = match(date_pattern, s)
    position_match = match(position_pattern, s)
    if isnothing(position_match)
        position_pattern = r"UTP \d+"
        position_match = match(position_pattern, s)
    end
    return date_match, position_match
end

function date_position_tag(s::String)
    date_match, position_match = parse_worm_dataset_path(s)
    if isnothing(date_match)
        if isnothing(position_match)
            error("Unable to parse date and position from $s")
        else
            return position_match.match
        end
    else
        if isnothing(position_match)
            return date_match.match
        else
            return string(date_match.match, "_", position_match.match)
        end
    end
    # not reachable
end

using HDF5

const annotations_cache = Dict{Tuple{String, UnitRange, Bool}, Vector}()
const annotation_position_cache = Dict{String, Any}()
#const my_annotation_position_cache = Dict{String, Any}()
const my_annotation_position_cache = Dict{String, Vector{Vector{Point3{Float64}}}}()

include("demo_averaging/load_straightened_annotations_over_time.jl")
include("demo_averaging/get_cell_trajectory_dict.jl")

include("demo_averaging/save_cache.jl")

# initialize my_annotation_position_cache
@info "Loading straightened annotation positions..."
load_annotation_cache()
@info "Loading warped annotation positions..."
load_annotations_cache()


function save_annotation_position_cache(
    filename::String,
    datasets::Dict{String, Vector{ShroffCelegansModels.Datasets.NormalizedDataset}},
    cache::Dict{String, Vector{Vector{Point3{Float64}}}};
    num_timepoints::Union{Symbol, Int} = :raw,
    expand_annotations::Bool = false,
    append_seam_cells::Bool = true
)
    h5open(filename, "w") do h5f
        for (group_name, embryos) in datasets
            h5g = create_group(h5f, group_name)
            println("group_name: $group_name")
            for embryo in embryos
                println("embryo.path: $(embryo.path)")
                try
                    tag = date_position_tag(embryo.path)
                    while haskey(h5g, tag)
                        error("Tag $tag already exists in $group_name")
                    end
                    # h5g[tag] = stack(cache[embryo.path])
                    # h5g[tag] = stack(get_cell_trajectory_dict(embryo; use_myuntwist = true))
                    #data = stack(get_cell_trajectory_dict(embryo; use_myuntwist = true))
                    dict = get_cell_trajectory_dict(embryo; use_myuntwist = true)
                    dict = filter(dict) do p
                        haskey(embryo.cell_key.mapping, Symbol(first(p)))
                    end
                    smts = StraightenedModelTimeSeries(embryo)
                    annotation_names = map(collect(keys(dict))) do name
                        # get(embryo.cell_key.mapping, Symbol(name), name)
                        embryo.cell_key.mapping[Symbol(name)]
                    end
                    if num_timepoints == :raw
                        num_timepoints = length(range(embryo.cell_key))
                    elseif num_timepoints isa Symbol
                        error("Invalid num_timepoints: $num_timepoints")
                    end
                    timepoints = LinRange(0, 1, num_timepoints)
                    data = map(values(dict)) do spline
                        if ismissing(spline)
                            fill(Point3{Float64}(NaN, NaN, NaN), num_timepoints)
                        else
                            spline.(timepoints)
                        end
                    end
                    if !isempty(data)
                        data = stack(data)
                        data = reshape(reinterpret(Float64, data), (3, size(data)...))
                    else
                        data = zeros(Float64, 3, num_timepoints, 0)
                    end
                    #data = reinterpret(Float64, data)
                    #data = reshape(data, (3, size(data)...))
                    # h5d = create_dataset(h5g, tag, datatype(Point3{Float64}), size(data))

                    # Add seam cells
                    if append_seam_cells
                        left_seam_cell_names = psuedo_seam_cells
                        right_seam_cell_names = replace.(left_seam_cell_names, 'L' => 'R')
                        seam_cell_names = vec(stack([left_seam_cell_names, right_seam_cell_names]; dims=1))
                        smts_length = length(range(embryo.cell_key))
                        seam_cell_data = map(timepoints) do ntp
                            # convert normalized time to timepoint
                            tp = ntp * (smts_length - 1) + 1
                            # straightened model
                            model = smts(tp)
                            if isnothing(model)
                                @info "Model is nothing" ntp tp smts_length
                            end
                            right_spline = transverse_spline(model, 1)
                            left_spline = transverse_spline(model, 17)
                            pts = BSplineKit.SplineInterpolations.interpolation_points(right_spline)
                            dict = Dict(Types.names(model) .=> vec(stack((left_spline.(pts), right_spline.(pts)); dims=1)))
                            seam_cell_locations = map(seam_cell_names) do name
                                if haskey(dict, name)
                                    dict[name]
                                else
                                    @warn "Seam cell $name not found $(embryo.path) $ntp $tp"
                                    Point3{Float64}(NaN, NaN, NaN)
                                end
                            end
                            return seam_cell_locations
                        end
                        seam_cell_data = stack(seam_cell_data; dims=1)
                        seam_cell_data = reshape(reinterpret(Float64, seam_cell_data), (3, size(seam_cell_data)...))
                        append!(annotation_names, seam_cell_names)
                        data = cat(data, seam_cell_data; dims=3)
                    end

                    if expand_annotations
                        h5g2 = create_group(h5g, tag)
                        for (i, name) in enumerate(annotation_names)
                            annotation_data = @view data[:,:,i]
                            h5d = create_dataset(h5g2, name, Float64, size(annotation_data))
                            write_dataset(h5d, datatype(Float64), annotation_data)
                        end
                        attributed_obj = h5g2
                    else
                        h5d = create_dataset(h5g, tag, Float64, size(data))
                        write_dataset(h5d, datatype(Float64), data)
                        attrs(h5d)["annotation_names"] = annotation_names
                        attributed_obj = h5d
                    end
                    attrs(attributed_obj)["path"] = embryo.path
                    attrs(attributed_obj)["normalized_time"] = collect(timepoints)
                    attrs(attributed_obj)["standard_minutes_post_fertilization"] = collect(timepoints .* 420 .+ 420)
                catch err
                    #display(err)
                    # println(err)
                    if err isa KeyError
                        println(err)
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end
end

function save_annotation_position_cache_all_dated(datasets::Dict{String, Vector{Datasets.NormalizedDataset}}; clear = true)
    # Clear caches
    if clear
        empty!(annotation_position_cache)
        empty!(my_annotation_position_cache)
        empty!(annotations_cache)
    end

    date_str = "$(Dates.today())"
    date_str = replace(date_str, "-" => "_")
    for (timepoints, expanded) in Iterators.product((:raw, 420), (true, false))
        tp_str = string(timepoints)
        if expanded
            tp_str *= "_expanded"
        end
        cache_file = "embryos_$(tp_str)_$date_str.h5"
        save_annotation_position_cache(
            cache_file,
            datasets,
            my_annotation_position_cache;
            num_timepoints = timepoints,
            expand_annotations = expanded,
            append_seam_cells = true
        )
        @info "Saved annotation position cache to $cache_file" timepoints expanded
    end
end