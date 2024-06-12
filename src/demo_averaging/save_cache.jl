function save_annotation_cache()
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

function load_annotation_cache()
    function _descend(p::Union{HDF5.File, HDF5.Group})
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

        _paths[1] = _paths[1] * ":\\\\"
        _path = joinpath(_paths...)

        data = d[]
    end
    h5open("my_annotation_position_cache.h5", "r") do h5f
        _descend(h5f)
    end
end