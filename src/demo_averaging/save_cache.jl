using HDF5
using Printf
using GeometryBasics

if !@isdefined(my_annotation_position_cache)
    const my_annotation_position_cache = Dict{String, Vector{Vector{Point3{Float64}}}}()
end
if !@isdefined(annotations_cache)
    const annotations_cache = Dict{Tuple{String, UnitRange, Bool}, Vector}()
end

function save_annotation_cache()
    # my_annotation_position_cache
    h5open("my_annotation_position_cache.h5", "w") do h5f
        for (k,v) in my_annotation_position_cache
            parts = splitpath(k)
            parts[1] = replace(parts[1], ":" => "", "\\" => "")
            group_name = join(parts, "/")
            for (idx, points) in pairs(v)
                _points = reinterpret(Float64, points)
                _points = reshape(_points, 3, :)
                _points = transpose(_points)
                h5f[group_name * "/" * @sprintf("%03d", idx)] = collect(_points)
            end
        end
    end
end

function save_annotations_cache(annotations_cache = annotations_cache)
    # annotations_cache
    h5open("annotations_cache.h5", "w") do h5f
        for (k,v) in annotations_cache
            _path, _range, _my_untwist = k
            parts = splitpath(_path)
            parts[1] = replace(parts[1], ":" => "", "\\" => "")
            group_name = join(parts, "/")
            h5g = create_group(h5f, group_name)
            attrs(h5g)["range_start"] = first(_range)
            attrs(h5g)["range_end"] = last(_range)
            attrs(h5g)["new_untwist"] = UInt8(_my_untwist)
            for (idx, data) in pairs(v)
                idx_str = @sprintf("%03d", idx)
                h5g_data = create_group(h5g, idx_str)
                if !ismissing(data)
                    for (k2, v2) in data
                        if k2 isa Integer
                            k2 = string(k2)
                        end
                        # Convert Point3{Float64} to Vector{Float64}
                        write_dataset(h5g_data, k2, collect(v2))
                        #h5g_data[k2] = collect(v2)
                    end
                end
            end
        end
    end
end

function load_annotations_cache(annotations_cache = annotations_cache)
    function _descend(p::Union{HDF5.File,HDF5.Group})
        for k in keys(p)
            _descend(p[k])
        end
    end
    function _descend(d::HDF5.Dataset)
        _name = HDF5.name(d)
        #println(_name)
        _paths = splitpath(_name)
        popfirst!(_paths)

        k2 = pop!(_paths)

        last_path = pop!(_paths)
        idx = tryparse(Int, last_path)
        P = nothing
        try
            P = parent(parent(d))
        catch err
            println(d)
            throw(err)
        end

        while isnothing(idx)
            k2 = last_path * "/" * k2
            last_path = pop!(_paths)
            idx = tryparse(Int, last_path)
            P = parent(P)
        end

        _paths[1] = _paths[1] * ":\\"
        _path = joinpath(_paths...)

        data = d[]
        pt = Point3{Float64}(data)
        #println(_path)
        #println(HDF5.name(d))
        _range_start = attrs(P)["range_start"]
        _range_end = attrs(P)["range_end"]
        _range = _range_start:_range_end

        cache = get!(annotations_cache, (_path, _range, true)) do
            N = length(_range)
            Vector{Union{Missing, Dict{String, Point3{Float64}}}}(missing, N)
        end
        data_cache = get(cache, idx, missing)
        if ismissing(data_cache)
            data_cache = Dict{String, Point3{Float64}}()
        end
        data_cache[k2] = pt
        cache[idx] = data_cache
    end
    h5open("annotations_cache.h5", "r") do h5f
        _descend(h5f)
    end

end

function load_annotation_cache()
    function _descend(p::Union{HDF5.File,HDF5.Group})
        for k in keys(p)
            _descend(p[k])
        end
    end
    function _descend(d::HDF5.Dataset)
        _name = HDF5.name(d)
        _paths = splitpath(_name)
        popfirst!(_paths)

        idx = pop!(_paths)
        idx = parse(Int, idx)

        _paths[1] = _paths[1] * ":\\"
        _path = joinpath(_paths...)

        data = d[]
        pts = Point3{Float64}.(eachrow(data))
        #println(_path)
        #println(HDF5.name(d))
        cache = get!(my_annotation_position_cache, _path) do
            P = parent(d)
            N = count(keys(P)) do k
                isa(P[k], HDF5.Dataset)
            end
            Vector{Vector{Point3{Float64}}}(undef, N)
        end
        cache[idx] = pts
    end
    h5open("my_annotation_position_cache.h5", "r") do h5f
        _descend(h5f)
    end
end
