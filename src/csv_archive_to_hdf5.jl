using HDF5: HDF5, h5open, attrs, create_group
using HDF5.Filters: Shuffle
using H5Zzstd: ZstdFilter
using CSV: CSV
using DataFrames: DataFrame, ncol, nrow
import DataFrames
using JSON3: JSON3
using Dates: Dates

const _CSV_ARCHIVE_SCHEMA_VERSION = 2
const _CSV_ARCHIVE_ZSTD_LEVEL = 9
# HDF5 chunked storage carries ~2 KB of per-dataset overhead. Skip the filter
# for arrays smaller than this.
const _CSV_ARCHIVE_MIN_COMPRESS_LEN = 64

const _POINT_SET_SPECS = (
    (name = "lattice",
     path = "lattice_final/lattice.csv",                              with_seg = false),
    (name = "annotations",
     path = "integrated_annotation/annotations.csv",                  with_seg = false),
    (name = "seam_cells",
     path = "seam_cell_final/seam_cells.csv",                         with_seg = false),
    (name = "named_seam_cells",
     path = "named_seam_cells/seam_cells.csv",                        with_seg = true),
    (name = "straightened_lattice",
     path = "straightened_lattice/straightened_lattice.csv",          with_seg = false),
    (name = "straightened_annotations",
     path = "straightened_annotations/straightened_annotations.csv",  with_seg = false),
    (name = "straightened_seam_cells",
     path = "straightened_seamcells/straightened_seamcells.csv",      with_seg = false),
    (name = "straightened_named_seam_cells",
     path = "straightened_named_seamcells/straightened_seamcells.csv", with_seg = false),
)

# ─── public entry point ──────────────────────────────────────────────────────

"""
    csv_archive_to_hdf5(archive_path, h5_path; overwrite=false)

Convert a per-dataset CSV tar.gz archive (as produced by
`web/scripts/gather_csv_archives.jl`) into a single Zstd-compressed HDF5 file.

To keep the file compact, **per-timepoint arrays are stacked into one tensor
per point-set type** (lattice, annotations, …) rather than producing a separate
group per timepoint. This collapses thousands of tiny datasets into a few
large compressible ones.

Layout (schema_version = $(_CSV_ARCHIVE_SCHEMA_VERSION)):

```
/                            attrs: archive_id, source_archive, created_at,
                                    schema_version, num_timepoints,
                                    dataset_name, dataset_start, dataset_end
├── timepoint_index          Int[T]    — the Decon_reg_<N> numbers, in order
├── cell_key/                attrs: end
│   ├── keys                 String[K]
│   └── values               String[K]
├── missed_segs_nuc/         (optional) names / volumes
├── missed_segs_seam/        (same)
├── point_sets/
│   └── <pset>/              one group per type in _POINT_SET_SPECS
│       ├── names            String[max_n] — union of names, first-encounter order
│       ├── counts           Int[T]        — active row count per timepoint
│       ├── xyz              Float64 (T, max_n, 3) NaN-padded
│       ├── rgb              UInt8   (T, max_n, 3) — when present
│       └── lattice_segment  Int     (T, max_n)    — named_seam_cells only
├── cross_sections/
│   ├── attrs: indices       Int[]   — union of indices seen
│   └── <i>/                 one group per cross-section index
│       ├── counts           Int[T]  — 0 when section i absent for that tp
│       └── xyz              Float64 (T, max_m, 3) NaN-padded
└── statistics/
    ├── diameters/           counts Int[T], data Float64 (T, max_D)
    ├── sample_planes/       attrs: columns; counts; data Float64 (T, max_S, 12)
    ├── frame_straight/      counts Int[T];
    │                        {center,left,right}_labels String[max_FS] (shared);
    │                        {center,left,right}_xyz Float64 (T, max_FS, 3)
    ├── lattice_info_before/ attrs: columns;
    │                        total_length Float64[T];
    │                        counts Int[T];
    │                        data Float64 (T, max_P, 4)
    └── lattice_info_after/  (same shape)
```

For each point-set tensor, slot `i` corresponds to the cell whose name is at
`names[i]`. If a timepoint doesn't include that cell, the slot is NaN-padded
(for `xyz`) or zero-filled (for `rgb` / `lattice_segment`). To find which
cells are present at timepoint `t`, use `!isnan.(xyz[t, :, 1])` — do *not*
slice by `counts[t]`, since active slots may be interleaved across the
union-ordered name list. `counts[t]` records the number of rows in the
original per-timepoint CSV (useful as a sanity check).

For `cross_sections/<i>/xyz` and `statistics/*/data`, rows ARE positional —
slice the first `counts[t]` rows to get a timepoint's active data.

Missing files for a given timepoint are silently skipped (`counts[t] = 0`).
"""
function csv_archive_to_hdf5(
    archive_path::AbstractString,
    h5_path::AbstractString;
    overwrite::Bool = false,
)
    isfile(archive_path) || error("Archive not found: $archive_path")
    if isfile(h5_path) && !overwrite
        error("HDF5 output already exists: $h5_path (pass overwrite=true to replace)")
    end

    mktempdir() do extract_dir
        run(`tar xzf $archive_path -C $extract_dir`)
        parent_dir, _, tp_dirs = _locate_dataset_dirs(extract_dir)

        snapshot = _parse_all_timepoints(tp_dirs)

        tmp_h5 = string(h5_path, ".tmp.", getpid(), ".", time_ns())
        try
            h5open(tmp_h5, "w") do h5
                attrs(h5)["archive_id"] = _archive_id(archive_path)
                attrs(h5)["source_archive"] = basename(archive_path)
                attrs(h5)["created_at"] = string(Dates.now(Dates.UTC))
                attrs(h5)["schema_version"] = _CSV_ARCHIVE_SCHEMA_VERSION
                attrs(h5)["num_timepoints"] = length(tp_dirs)

                _write_parent_files!(h5, parent_dir)
                _write_zstd!(h5, "timepoint_index", snapshot.tp_numbers)

                _write_stacked_point_sets!(h5, snapshot)
                _write_stacked_cross_sections!(h5, snapshot)
                _write_stacked_statistics!(h5, snapshot)
            end
            mv(tmp_h5, h5_path; force = true)
        catch
            isfile(tmp_h5) && rm(tmp_h5; force = true)
            rethrow()
        end
    end
    return h5_path
end

# ─── write helper ────────────────────────────────────────────────────────────

"""
    _write_zstd!(group, name, data)

Write `data` to `group[name]`. Numeric arrays of length ≥
`_CSV_ARCHIVE_MIN_COMPRESS_LEN` get chunked storage with Shuffle + Zstd.
Strings and tiny arrays are written inline (no filter).
"""
function _write_zstd!(g, name::AbstractString, data)
    if eltype(data) <: AbstractString || length(data) < _CSV_ARCHIVE_MIN_COMPRESS_LEN
        g[name] = data
        return nothing
    end
    sz = size(data)
    filters = if eltype(data) <: Union{AbstractFloat,Integer}
        [Shuffle(), ZstdFilter(_CSV_ARCHIVE_ZSTD_LEVEL)]
    else
        [ZstdFilter(_CSV_ARCHIVE_ZSTD_LEVEL)]
    end
    d = HDF5.create_dataset(
        g, name, eltype(data), sz;
        chunk = sz, filters = filters,
    )
    d[fill(:, length(sz))...] = data
    return nothing
end

# ─── dataset-root discovery ──────────────────────────────────────────────────

function _locate_dataset_dirs(extract_dir::AbstractString)
    tp_dirs = String[]
    for (root, dirs, _) in walkdir(extract_dir)
        for d in dirs
            occursin(r"^Decon_reg_\d+_results$", d) && push!(tp_dirs, joinpath(root, d))
        end
    end
    isempty(tp_dirs) && error("Archive contains no Decon_reg_*_results timepoint folders")

    regb_candidates = unique(dirname(dirname(d)) for d in tp_dirs)
    length(regb_candidates) == 1 ||
        error("Ambiguous dataset layout — multiple RegB-equivalent dirs: $regb_candidates")
    regb_dir = first(regb_candidates)
    parent_dir = dirname(regb_dir)

    sort!(tp_dirs; by = _timepoint_number)
    return parent_dir, regb_dir, tp_dirs
end

function _timepoint_number(tp_dir::AbstractString)
    m = match(r"Decon_reg_(\d+)_results$", basename(tp_dir))
    m === nothing && error("Unparseable timepoint folder name: $tp_dir")
    return parse(Int, m.captures[1])
end

function _archive_id(archive_path::AbstractString)
    b = basename(archive_path)
    endswith(b, ".tar.gz") && return b[1:(end - length(".tar.gz"))]
    endswith(b, ".tgz") && return b[1:(end - length(".tgz"))]
    return b
end

# ─── parent-level files (CellKey, cell_key.json, MissedSegs_*) ───────────────

function _write_parent_files!(h5, parent_dir::AbstractString)
    cellkey_csv = joinpath(parent_dir, "CellKey.csv")
    if isfile(cellkey_csv)
        name, start_tp, end_tp = _parse_cellkey_csv(cellkey_csv)
        name === nothing || (attrs(h5)["dataset_name"] = name)
        start_tp === nothing || (attrs(h5)["dataset_start"] = start_tp)
        end_tp === nothing || (attrs(h5)["dataset_end"] = end_tp)
    end

    cellkey_json = joinpath(parent_dir, "cell_key.json")
    if isfile(cellkey_json)
        ks, vs, json_end = _parse_cell_key_json(cellkey_json)
        ck = create_group(h5, "cell_key")
        _write_zstd!(ck, "keys", ks)
        _write_zstd!(ck, "values", vs)
        json_end === nothing || (attrs(ck)["end"] = json_end)
    end

    _write_missed_segs!(h5, joinpath(parent_dir, "MissedSegs_nuc.csv"), "missed_segs_nuc")
    _write_missed_segs!(h5, joinpath(parent_dir, "MissedSegs_seam.csv"), "missed_segs_seam")
    return nothing
end

function _parse_cellkey_csv(path::AbstractString)
    name = nothing; start_tp = nothing; end_tp = nothing
    try
        lines = readlines(path)
        if length(lines) >= 1
            parts = split(lines[1], ',')
            !isempty(parts) && !isempty(strip(parts[1])) && (name = String(strip(parts[1])))
        end
        if length(lines) >= 2
            parts = split(lines[2], ',')
            if length(parts) >= 2
                s = tryparse(Int, strip(parts[1]))
                e = tryparse(Int, strip(parts[2]))
                s === nothing || (start_tp = s)
                e === nothing || (end_tp = e)
            end
        end
    catch err
        @warn "Failed to parse CellKey.csv" path err
    end
    return name, start_tp, end_tp
end

function _parse_cell_key_json(path::AbstractString)
    obj = JSON3.read(read(path, String))
    mapping = haskey(obj, :mapping) ? obj.mapping : obj
    ks = String[]; vs = String[]
    for (k, v) in pairs(mapping)
        push!(ks, String(k)); push!(vs, String(v))
    end
    json_end = haskey(obj, :end) ? Int(obj.end) : nothing
    return ks, vs, json_end
end

function _write_missed_segs!(h5, csv_path::AbstractString, group_name::AbstractString)
    isfile(csv_path) || return nothing
    df = try
        CSV.read(csv_path, DataFrame; header = 2)
    catch err
        @warn "Failed to parse MissedSegs file" csv_path err
        return nothing
    end
    ncol(df) >= 2 || return nothing
    g = create_group(h5, group_name)
    _write_zstd!(g, "names", String.(string.(df[!, 1])))
    _write_zstd!(g, "volumes", _as_float64(df[!, 2]))
    return nothing
end

# ─── timepoint parser (in-memory snapshot) ───────────────────────────────────

struct TimepointPointSet
    names::Vector{String}
    xyz::Matrix{Float64}                  # (n, 3)
    rgb::Union{Nothing,Matrix{UInt8}}     # (n, 3) or nothing
    lattice_segment::Union{Nothing,Vector{Int}}
end

struct TimepointFrameStraight
    center_labels::Vector{String}
    left_labels::Vector{String}
    right_labels::Vector{String}
    center_xyz::Matrix{Float64}           # (n, 3)
    left_xyz::Matrix{Float64}             # (n, 3)
    right_xyz::Matrix{Float64}            # (n, 3)
end

struct TimepointLatticeInfo
    total_length::Union{Nothing,Float64}
    data::Matrix{Float64}                 # (p, 4)
end

struct TimepointData
    tp_number::Int
    point_sets::Dict{String,TimepointPointSet}
    cross_sections::Dict{Int,Matrix{Float64}}   # idx → (m, 3)
    diameters::Union{Nothing,Vector{Float64}}
    sample_planes::Union{Nothing,Matrix{Float64}}   # (s, 12)
    frame_straight::Union{Nothing,TimepointFrameStraight}
    lattice_info_before::Union{Nothing,TimepointLatticeInfo}
    lattice_info_after::Union{Nothing,TimepointLatticeInfo}
end

struct Snapshot
    tp_numbers::Vector{Int}
    timepoints::Vector{TimepointData}
end

function _parse_all_timepoints(tp_dirs::Vector{String})
    tps = TimepointData[]
    for tp_dir in tp_dirs
        push!(tps, _parse_one_timepoint(tp_dir))
    end
    return Snapshot([t.tp_number for t in tps], tps)
end

function _parse_one_timepoint(tp_dir::AbstractString)
    point_sets = Dict{String,TimepointPointSet}()
    for spec in _POINT_SET_SPECS
        path = joinpath(tp_dir, spec.path)
        isfile(path) || continue
        ps = _parse_point_set(path; with_lattice_segment = spec.with_seg)
        ps === nothing || (point_sets[spec.name] = ps)
    end

    cross_sections = Dict{Int,Matrix{Float64}}()
    cs_dir = joinpath(tp_dir, "model_crossSections")
    if isdir(cs_dir)
        for f in readdir(cs_dir)
            m = match(r"^latticeCrossSection_(\d+)\.csv$", f)
            m === nothing && continue
            idx = parse(Int, m.captures[1])
            mat = _parse_cross_section(joinpath(cs_dir, f))
            mat === nothing || (cross_sections[idx] = mat)
        end
    end

    diameters = nothing
    sample_planes = nothing
    frame_straight = nothing
    lattice_info_before = nothing
    lattice_info_after = nothing
    stats_dir = joinpath(tp_dir, "statistics")
    if isdir(stats_dir)
        diameters = _parse_diameters(joinpath(stats_dir, "Diameters.csv"))
        sample_planes = _parse_sample_planes(joinpath(stats_dir, "SamplePlanes.csv"))
        fs = _find_frame_straight(stats_dir)
        fs === nothing || (frame_straight = _parse_frame_straight(fs))
        lattice_info_before = _parse_lattice_info(joinpath(stats_dir, "LatticeInfo_before.csv"))
        lattice_info_after = _parse_lattice_info(joinpath(stats_dir, "LatticeInfo_after.csv"))
    end

    return TimepointData(
        _timepoint_number(tp_dir),
        point_sets,
        cross_sections,
        diameters,
        sample_planes,
        frame_straight,
        lattice_info_before,
        lattice_info_after,
    )
end

function _parse_point_set(csv_path::AbstractString; with_lattice_segment::Bool = false)
    df = try
        CSV.read(csv_path, DataFrame)
    catch err
        @warn "Failed to parse point-set CSV" csv_path err
        return nothing
    end
    name_col = _find_column(df, ["name"])
    x_col = _find_column(df, ["x_voxels", "x"])
    y_col = _find_column(df, ["y_voxels", "y"])
    z_col = _find_column(df, ["z_voxels", "z"])
    (name_col === nothing || x_col === nothing || y_col === nothing || z_col === nothing) &&
        return nothing

    names_v = String.(string.(df[!, name_col]))
    xyz = _xyz_matrix(df, x_col, y_col, z_col)
    rgb = nothing
    r = _find_column(df, ["R"]); gc = _find_column(df, ["G"]); b = _find_column(df, ["B"])
    if r !== nothing && gc !== nothing && b !== nothing
        rgb = hcat(_as_uint8(df[!, r]), _as_uint8(df[!, gc]), _as_uint8(df[!, b]))
    end
    seg = nothing
    if with_lattice_segment
        seg_col = _find_column(df, ["lattice segment", "lattice_segment"])
        seg_col === nothing || (seg = _as_int(df[!, seg_col]))
    end
    return TimepointPointSet(names_v, xyz, rgb, seg)
end

function _parse_cross_section(csv_path::AbstractString)
    df = try
        CSV.read(csv_path, DataFrame; header = 2)
    catch err
        @warn "Failed to parse crossSection CSV" csv_path err
        return nothing
    end
    ncol(df) >= 3 || return nothing
    return hcat(_as_float64(df[!, 1]), _as_float64(df[!, 2]), _as_float64(df[!, 3]))
end

function _parse_diameters(csv_path::AbstractString)
    isfile(csv_path) || return nothing
    df = try
        CSV.read(csv_path, DataFrame)
    catch err
        @warn "Failed to parse Diameters.csv" csv_path err
        return nothing
    end
    ncol(df) >= 1 || return nothing
    return _as_float64(df[!, 1])
end

function _parse_sample_planes(csv_path::AbstractString)
    isfile(csv_path) || return nothing
    df = try
        CSV.read(csv_path, DataFrame; header = false, skipto = 2)
    catch err
        @warn "Failed to parse SamplePlanes.csv" csv_path err
        return nothing
    end
    ncol(df) == 12 || return nothing
    mat = Matrix{Float64}(undef, nrow(df), 12)
    for j in 1:12
        mat[:, j] = _as_float64(df[!, j])
    end
    return mat
end

function _find_frame_straight(stats_dir::AbstractString)
    for f in readdir(stats_dir)
        endswith(f, "_Frame_Straight.csv") && return joinpath(stats_dir, f)
    end
    return nothing
end

function _parse_frame_straight(csv_path::AbstractString)
    df = try
        CSV.read(csv_path, DataFrame; header = false, skipto = 2)
    catch err
        @warn "Failed to parse Frame_Straight.csv" csv_path err
        return nothing
    end
    ncol(df) == 12 || return nothing
    return TimepointFrameStraight(
        String.(string.(df[!, 1])),
        String.(string.(df[!, 5])),
        String.(string.(df[!, 9])),
        hcat(_as_float64(df[!, 2]),  _as_float64(df[!, 3]),  _as_float64(df[!, 4])),
        hcat(_as_float64(df[!, 6]),  _as_float64(df[!, 7]),  _as_float64(df[!, 8])),
        hcat(_as_float64(df[!, 10]), _as_float64(df[!, 11]), _as_float64(df[!, 12])),
    )
end

function _parse_lattice_info(csv_path::AbstractString)
    isfile(csv_path) || return nothing
    lines = readlines(csv_path)
    total_length = nothing
    if !isempty(lines)
        parts = split(lines[1], ',')
        if length(parts) >= 2
            v = tryparse(Float64, strip(parts[2]))
            v === nothing || (total_length = v)
        end
    end
    header_idx = findfirst(l -> startswith(strip(l), "pair"), lines)
    header_idx === nothing && return nothing
    df = try
        CSV.read(IOBuffer(join(lines[header_idx:end], "\n")), DataFrame)
    catch err
        @warn "Failed to parse LatticeInfo CSV" csv_path err
        return nothing
    end
    ncol(df) >= 4 || return nothing
    mat = Matrix{Float64}(undef, nrow(df), 4)
    for j in 1:4
        mat[:, j] = _as_float64(df[!, j])
    end
    return TimepointLatticeInfo(total_length, mat)
end

# ─── stacked writers (the heart of the size reduction) ───────────────────────

function _write_stacked_point_sets!(h5, snap::Snapshot)
    isempty(snap.timepoints) && return nothing
    ps_root = create_group(h5, "point_sets")
    T = length(snap.timepoints)

    for spec in _POINT_SET_SPECS
        # Skip if no timepoint produced this point set
        any(haskey(tp.point_sets, spec.name) for tp in snap.timepoints) || continue

        # Build union of names in first-encounter order.
        seen = String[]
        seen_set = Set{String}()
        for tp in snap.timepoints
            ps = get(tp.point_sets, spec.name, nothing)
            ps === nothing && continue
            for n in ps.names
                if !(n in seen_set)
                    push!(seen, n); push!(seen_set, n)
                end
            end
        end
        max_n = length(seen)
        max_n == 0 && continue
        name_to_idx = Dict(seen[i] => i for i in 1:max_n)

        has_rgb = any(
            haskey(tp.point_sets, spec.name) &&
            tp.point_sets[spec.name].rgb !== nothing
            for tp in snap.timepoints
        )
        has_seg = spec.with_seg && any(
            haskey(tp.point_sets, spec.name) &&
            tp.point_sets[spec.name].lattice_segment !== nothing
            for tp in snap.timepoints
        )

        xyz = fill(NaN, T, max_n, 3)
        rgb = has_rgb ? zeros(UInt8, T, max_n, 3) : nothing
        seg = has_seg ? zeros(Int, T, max_n) : nothing
        counts = zeros(Int, T)

        for (t, tp) in enumerate(snap.timepoints)
            ps = get(tp.point_sets, spec.name, nothing)
            ps === nothing && continue
            counts[t] = length(ps.names)
            for (row_in_tp, name) in enumerate(ps.names)
                slot = name_to_idx[name]
                @views xyz[t, slot, :] .= ps.xyz[row_in_tp, :]
                if has_rgb && ps.rgb !== nothing
                    @views rgb[t, slot, :] .= ps.rgb[row_in_tp, :]
                end
                if has_seg && ps.lattice_segment !== nothing
                    seg[t, slot] = ps.lattice_segment[row_in_tp]
                end
            end
        end

        g = create_group(ps_root, spec.name)
        _write_zstd!(g, "names", seen)
        _write_zstd!(g, "counts", counts)
        _write_zstd!(g, "xyz", xyz)
        has_rgb && _write_zstd!(g, "rgb", rgb)
        has_seg && _write_zstd!(g, "lattice_segment", seg)
    end
    return nothing
end

function _write_stacked_cross_sections!(h5, snap::Snapshot)
    T = length(snap.timepoints)
    all_indices = sort!(collect(reduce(union,
        (keys(tp.cross_sections) for tp in snap.timepoints); init = Set{Int}())))
    isempty(all_indices) && return nothing

    cs_root = create_group(h5, "cross_sections")
    attrs(cs_root)["indices"] = all_indices

    for idx in all_indices
        max_m = 0
        for tp in snap.timepoints
            haskey(tp.cross_sections, idx) || continue
            max_m = max(max_m, size(tp.cross_sections[idx], 1))
        end
        xyz = fill(NaN, T, max_m, 3)
        counts = zeros(Int, T)
        for (t, tp) in enumerate(snap.timepoints)
            haskey(tp.cross_sections, idx) || continue
            mat = tp.cross_sections[idx]
            m = size(mat, 1)
            counts[t] = m
            @views xyz[t, 1:m, :] .= mat
        end
        g = create_group(cs_root, string(idx))
        _write_zstd!(g, "counts", counts)
        _write_zstd!(g, "xyz", xyz)
    end
    return nothing
end

function _write_stacked_statistics!(h5, snap::Snapshot)
    T = length(snap.timepoints)
    stats = create_group(h5, "statistics")

    # diameters: (T, max_D) Float64, counts
    if any(tp.diameters !== nothing for tp in snap.timepoints)
        max_d = maximum(
            tp.diameters === nothing ? 0 : length(tp.diameters) for tp in snap.timepoints
        )
        data = fill(NaN, T, max_d)
        counts = zeros(Int, T)
        for (t, tp) in enumerate(snap.timepoints)
            tp.diameters === nothing && continue
            n = length(tp.diameters)
            counts[t] = n
            @views data[t, 1:n] .= tp.diameters
        end
        g = create_group(stats, "diameters")
        _write_zstd!(g, "counts", counts)
        _write_zstd!(g, "data", data)
    end

    # sample_planes: (T, max_S, 12) Float64
    if any(tp.sample_planes !== nothing for tp in snap.timepoints)
        max_s = maximum(
            tp.sample_planes === nothing ? 0 : size(tp.sample_planes, 1)
            for tp in snap.timepoints
        )
        data = fill(NaN, T, max_s, 12)
        counts = zeros(Int, T)
        for (t, tp) in enumerate(snap.timepoints)
            tp.sample_planes === nothing && continue
            n = size(tp.sample_planes, 1)
            counts[t] = n
            @views data[t, 1:n, :] .= tp.sample_planes
        end
        g = create_group(stats, "sample_planes")
        attrs(g)["columns"] =
            ["X1","Y1","Z1","X2","Y2","Z2","X3","Y3","Z3","X4","Y4","Z4"]
        _write_zstd!(g, "counts", counts)
        _write_zstd!(g, "data", data)
    end

    # frame_straight: stacked center/left/right xyz, shared labels
    if any(tp.frame_straight !== nothing for tp in snap.timepoints)
        max_fs = maximum(
            tp.frame_straight === nothing ? 0 : size(tp.frame_straight.center_xyz, 1)
            for tp in snap.timepoints
        )
        # Pick shared labels from the first timepoint that has them.
        ref_fs = nothing
        for tp in snap.timepoints
            tp.frame_straight === nothing && continue
            if size(tp.frame_straight.center_xyz, 1) == max_fs
                ref_fs = tp.frame_straight
                break
            end
            ref_fs === nothing && (ref_fs = tp.frame_straight)
        end

        c_xyz = fill(NaN, T, max_fs, 3)
        l_xyz = fill(NaN, T, max_fs, 3)
        r_xyz = fill(NaN, T, max_fs, 3)
        counts = zeros(Int, T)
        for (t, tp) in enumerate(snap.timepoints)
            tp.frame_straight === nothing && continue
            fs = tp.frame_straight
            n = size(fs.center_xyz, 1)
            counts[t] = n
            @views c_xyz[t, 1:n, :] .= fs.center_xyz
            @views l_xyz[t, 1:n, :] .= fs.left_xyz
            @views r_xyz[t, 1:n, :] .= fs.right_xyz
        end
        g = create_group(stats, "frame_straight")
        _write_zstd!(g, "counts", counts)
        if ref_fs !== nothing
            _write_zstd!(g, "center_labels", _pad_labels(ref_fs.center_labels, max_fs))
            _write_zstd!(g, "left_labels", _pad_labels(ref_fs.left_labels, max_fs))
            _write_zstd!(g, "right_labels", _pad_labels(ref_fs.right_labels, max_fs))
        end
        _write_zstd!(g, "center_xyz", c_xyz)
        _write_zstd!(g, "left_xyz", l_xyz)
        _write_zstd!(g, "right_xyz", r_xyz)
    end

    _write_stacked_lattice_info!(stats, snap, "lattice_info_before",
        tp -> tp.lattice_info_before)
    _write_stacked_lattice_info!(stats, snap, "lattice_info_after",
        tp -> tp.lattice_info_after)
    return nothing
end

function _write_stacked_lattice_info!(stats, snap::Snapshot, group_name::AbstractString,
                                       getter)
    T = length(snap.timepoints)
    any(getter(tp) !== nothing for tp in snap.timepoints) || return nothing

    max_p = maximum(
        getter(tp) === nothing ? 0 : size(getter(tp).data, 1) for tp in snap.timepoints
    )
    data = fill(NaN, T, max_p, 4)
    total_lengths = fill(NaN, T)
    counts = zeros(Int, T)
    for (t, tp) in enumerate(snap.timepoints)
        li = getter(tp)
        li === nothing && continue
        n = size(li.data, 1)
        counts[t] = n
        @views data[t, 1:n, :] .= li.data
        li.total_length === nothing || (total_lengths[t] = li.total_length)
    end
    g = create_group(stats, group_name)
    attrs(g)["columns"] = ["pair", "diameter", "left_distance", "right_distance"]
    _write_zstd!(g, "counts", counts)
    _write_zstd!(g, "total_length", total_lengths)
    _write_zstd!(g, "data", data)
    return nothing
end

# ─── column-coercion helpers ─────────────────────────────────────────────────

function _find_column(df::DataFrame, candidates)
    ns = string.(DataFrames.names(df))
    for c in candidates
        idx = findfirst(==(c), ns)
        idx === nothing || return idx
    end
    return nothing
end

function _xyz_matrix(df::DataFrame, x_col::Integer, y_col::Integer, z_col::Integer)
    n = size(df, 1)
    mat = Matrix{Float64}(undef, n, 3)
    mat[:, 1] = _as_float64(df[!, x_col])
    mat[:, 2] = _as_float64(df[!, y_col])
    mat[:, 3] = _as_float64(df[!, z_col])
    return mat
end

function _pad_labels(labels::Vector{String}, max_n::Integer)
    length(labels) >= max_n && return labels[1:max_n]
    out = copy(labels)
    while length(out) < max_n
        push!(out, "")
    end
    return out
end

_as_float64(v) = [Float64(x) for x in v]
_as_uint8(v) = [UInt8(clamp(Int(x), 0, 255)) for x in v]
_as_int(v) = [Int(x) for x in v]
