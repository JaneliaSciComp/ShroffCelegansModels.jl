using HDF5
using Printf
using GeometryBasics

if !@isdefined(my_annotation_position_cache)
    const my_annotation_position_cache = Dict{String, Vector{Vector{Point3{Float64}}}}()
end

struct AnnotationsCacheValue
    annotations::Vector{Union{Missing, Dict{String, Point3{Float64}}}}
    mtime::Float64
end

if !@isdefined(annotations_cache)
    const annotations_cache = Dict{Tuple{String, UnitRange, Bool}, AnnotationsCacheValue}()
end

function save_annotation_cache(; filename = "my_annotation_position_cache.h5")
    # my_annotation_position_cache
    h5open(filename, "w") do h5f
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

function save_annotations_cache(
    annotations_cache = annotations_cache;
    filename = joinpath(@__DIR__, "..", "..", "annotations_cache.h5")
)
    # annotations_cache
    h5open(filename, "w") do h5f
        for (k,v) in annotations_cache
            _path, _range, _my_untwist = k
            parts = splitpath(_path)
            parts[1] = replace(parts[1], ":" => "", "\\" => "")
            group_name = join(parts, "/")
            h5g = create_group(h5f, group_name)
            attrs(h5g)["range_start"] = first(_range)
            attrs(h5g)["range_end"] = last(_range)
            attrs(h5g)["new_untwist"] = UInt8(_my_untwist)
            per_idx_mtimes = _stat_per_idx_mtimes(_path, _range, length(v.annotations))
            for (idx, data) in pairs(v.annotations)
                idx_str = @sprintf("%03d", idx)
                h5g_data = create_group(h5g, idx_str)
                attrs(h5g_data)["mtime"] = per_idx_mtimes[idx]
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

# Stat integrated_annotation/annotations.csv per timepoint and return unix
# mtimes (NaN if the file doesn't exist or can't be stat'd). `_range` is the
# dataset's timepoint range, so the actual timepoint for offset i is _range[i].
function _stat_per_idx_mtimes(stored_path::AbstractString, _range::AbstractUnitRange{Int}, n::Int)::Vector{Float64}
    local_path = stored_path
    if Sys.isunix()
        # Tolerate either "X:\foo\bar" or "X:\\foo\\bar"; convert to nearline.
        local_path = replace(local_path, r"^[A-Za-z]:\\+" => "/nearline/shroff/")
        local_path = replace(local_path, "\\" => "/")
    end
    mtimes = fill(NaN, n)
    for i in 1:n
        timepoint = _range[i]
        try
            filepath = joinpath(
                local_path,
                "Decon_reg_$(timepoint)",
                "Decon_reg_$(timepoint)_results",
                "integrated_annotation",
                "annotations.csv",
            )
            if isfile(filepath)
                mtimes[i] = stat(filepath).mtime
            end
        catch
        end
    end
    return mtimes
end

function load_annotations_cache(
    annotations_cache = annotations_cache;
    filename = joinpath(@__DIR__, "..", "..", "annotations_cache.h5")
)
    if !isfile(filename)
        @warn "Annotations cache not found; leaving annotations_cache empty (will be computed lazily)" filename
        return annotations_cache
    end
    # Accumulate per-key state during traversal, then build immutable
    # AnnotationsCacheValue entries in a single finalize pass.
    annotations_by_key = Dict{
        Tuple{String, UnitRange{Int}, Bool},
        Vector{Union{Missing, Dict{String, Point3{Float64}}}}
    }()
    mtimes_by_key = Dict{Tuple{String, UnitRange{Int}, Bool}, Vector{Float64}}()

    function _descend(p::Union{HDF5.File,HDF5.Group})
        for k in keys(p)
            _descend(p[k])
        end
    end
    function _descend(d::HDF5.Dataset)
        _name = HDF5.name(d)
        _paths = split(_name, "/")
        popfirst!(_paths)

        k2 = pop!(_paths)
        last_path = pop!(_paths)
        idx = tryparse(Int, last_path)

        # idx_group is the timepoint group (where the mtime attr lives).
        # path_group is the dataset group (where range_start/_end live).
        idx_group = nothing
        path_group = nothing
        try
            idx_group = parent(d)
            path_group = parent(idx_group)
        catch err
            @warn "Could not load $d" err
            return
        end

        # When the annotation name contains '/', the leaf dataset is nested
        # inside extra subgroups; walk up until last_path parses as the idx.
        while isnothing(idx)
            k2 = last_path * "/" * k2
            last_path = pop!(_paths)
            idx = tryparse(Int, last_path)
            idx_group = path_group
            path_group = parent(path_group)
        end

        _paths[1] = _paths[1] * ":"
        _path = join(_paths, "\\")

        data = d[]
        pt = Point3{Float64}(data)

        _range_start = Int(attrs(path_group)["range_start"])
        _range_end = Int(attrs(path_group)["range_end"])
        _range = _range_start:_range_end
        key = (_path, _range, true)

        annotations = get!(annotations_by_key, key) do
            N = length(_range)
            Vector{Union{Missing, Dict{String, Point3{Float64}}}}(missing, N)
        end
        if ismissing(annotations[idx])
            annotations[idx] = Dict{String, Point3{Float64}}()
        end
        annotations[idx][k2] = pt

        mtimes = get!(mtimes_by_key, key) do
            fill(NaN, length(_range))
        end
        # Backward-compat: older files have no "mtime" attr → leave as NaN.
        try
            if haskey(attrs(idx_group), "mtime")
                mtimes[idx] = Float64(read(attrs(idx_group)["mtime"]))
            end
        catch err
            @warn "Could not read mtime attr for $(HDF5.name(idx_group))" err
        end
    end

    @info "Loading annotations cache from $filename"
    h5open(filename, "r") do h5f
        _descend(h5f)
    end

    # Reduce per-idx mtimes to a single dataset-level max (NaN if all missing).
    for (key, annotations) in annotations_by_key
        mtimes = mtimes_by_key[key]
        finite = filter(!isnan, mtimes)
        mtime = isempty(finite) ? NaN : maximum(finite)
        annotations_cache[key] = AnnotationsCacheValue(annotations, mtime)
    end

    return annotations_cache
end

function load_annotation_cache(; filename = joinpath(@__DIR__, "..", "..", "my_annotation_position_cache.h5"))
    if !isfile(filename)
        @warn "Annotation position cache not found; leaving my_annotation_position_cache empty (will be computed lazily)" filename
        return my_annotation_position_cache
    end
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

        # Reconstruct the cache key from the HDF5 group hierarchy. A single-
        # character root is a Windows drive letter (e.g. "X") written from a
        # Windows-built cache — keep the "X:\…" form so alias_cache_unix maps it
        # to the Linux dataset path. Any longer root (e.g. "nearline", written by
        # the recompute pipeline on Linux) is already an absolute Linux path that
        # equals dataset.path, so reconstruct it directly without mangling.
        if length(_paths[1]) == 1
            _paths[1] = _paths[1] * ":\\"
            _path = joinpath(_paths...)
        else
            _path = "/" * join(_paths, "/")
        end

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
    h5open(filename, "r") do h5f
        _descend(h5f)
    end
    return my_annotation_position_cache
end
