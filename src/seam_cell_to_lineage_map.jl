using ShroffCelegansModels: CSV, DataFrames, GeometryBasics, HDF5
using ShroffCelegansModels: transverse_splines
using CSV: CSV
using DataFrames: DataFrame, transform!, subset, select, nrow, eachrow, ByRow
using Printf: @sprintf
using LinearAlgebra: normalize, cross, dot, norm, I
#using GeometryBasics
#using GLMakie

# `save_celegans_avg_models.jl`, `seam_cell_pts.jl`, and `modelio.jl`
# are already included by the package (src/ShroffCelegansModels.jl and
# its transitive includes); no need to include them again here.

# Default search paths for the ryan_data CSVs the helpers below use.
# Resolved at @__DIR__ time so the package can be loaded from any cwd.
const _RYAN_DATA_DIR = joinpath(@__DIR__, "..", "ryan_data")
const _PRETWITCH_CSV = joinpath(_RYAN_DATA_DIR, "final_20251205_pre-twitch_coords.csv")
const _NAMING_CORRELATIONS_CSV = joinpath(_RYAN_DATA_DIR, "MIPAV_PositionalModel_Packer_Naming_Correlations_v6.csv")
const _COLOR_CODE_CSV = joinpath(_RYAN_DATA_DIR, "Updated_Color_Code_Assignments_04302026.csv")

# Mapping from embryonic lineage strings to seam cell names
const lineage_to_seam_cell_map = Dict(
    "ABplaaappa" => "H0L",
    "ABarpapppa" => "H0R",
    "ABplaaappp" => "H1L",
    "ABarpapppp" => "H1R",
    "ABarppaaap" => "H2L",
    "ABarpppaap" => "H2R",
    "ABarppapaa" => "V1L",
    "ABarppppaa" => "V1R",
    "ABarppapap" => "V2L",
    "ABarppppap" => "V2R",
    "ABplappapa" => "V3L",
    "ABprappapa" => "V3R",
    "ABarppappa" => "V4L",
    "ABarpppppa" => "V4R",
    "ABplapapaap" => "V5L",
    "ABprapapaap" => "V5R",
    "ABarppappp" => "V6L",
    "ABarpppppp" => "V6R",
    "ABplappppp" => "TL",
    "ABprappppp" => "TR"
)

# Reverse mapping from seam cell names to embryonic lineage strings
const seam_cell_to_lineage_map = Dict(
    "H0L" => "ABplaaappa",
    "H0R" => "ABarpapppa",
    "H1L" => "ABplaaappp",
    "H1R" => "ABarpapppp",
    "H2L" => "ABarppaaap",
    "H2R" => "ABarpppaap",
    "V1L" => "ABarppapaa",
    "V1R" => "ABarppppaa",
    "V2L" => "ABarppapap",
    "V2R" => "ABarppppap",
    "V3L" => "ABplappapa",
    "V3R" => "ABprappapa",
    "V4L" => "ABarppappa",
    "V4R" => "ABarpppppa",
    "V5L" => "ABplapapaap",
    "V5R" => "ABprapapaap",
    "V6L" => "ABarppappp",
    "V6R" => "ABarpppppp",
    "TL" => "ABplappppp",
    "TR" => "ABprappppp"
)

const left_seam_cells = ["H0L", "H1L", "H2L", "V1L", "V2L", "V3L", "V4L", "V5L", "V6L", "TL"]

const right_seam_cells = ["H0R", "H1R", "H2R", "V1R", "V2R", "V3R", "V4R", "V5R", "V6R", "TR"]

const seam_cell_midpoints = replace.(left_seam_cells, 'L' => 'M')

function get_pretwitch_df(path::AbstractString = _PRETWITCH_CSV)
    pretwitch_df = CSV.read(path, DataFrame)
    transform!(pretwitch_df, :cell => ByRow(strip) => :cell)
    return pretwitch_df
end

function get_lineage_df(lineage::Union{String,SubString}; pretwitch_df=get_pretwitch_df())
    return subset(pretwitch_df, :cell => cell -> strip.(cell) .== lineage)
end

function get_full_lineage_df(lineage::String; pretwitch_df=get_pretwitch_df())
    dfs = DataFrame[]
    for i in 3:length(lineage)
        sub_lineage = @view lineage[1:i]
        df = get_lineage_df(sub_lineage; pretwitch_df=pretwitch_df)
        push!(dfs, df)
    end
    return vcat(dfs...)
end

function get_seam_cell_df(seam_cell::String; pretwitch_df=get_pretwitch_df())
    lineage = seam_cell_to_lineage_map[seam_cell]
    return get_lineage_df(lineage; pretwitch_df=pretwitch_df)
end

function get_seam_cell_time_and_points(seam_cell::String)
    df = get_seam_cell_df(seam_cell)
    times = df.time .|> Float64
    points = [Point3(df.x[i], df.y[i], df.z[i]) for i in 1:nrow(df)]
    return times, points
end

function scatter_times_and_points!(times, points)
    scatter!(points, color=times)
end

function scatter_times_and_points(times, points)
    scatter(points, color=times)
end

function dataframe_to_point3f_array(df::DataFrame)
    return Point3f.(Vector.(eachrow(select(df, [:x, :y, :z]))))
end

function get_pretwitch_points_at_time(pretwitch_df, time::Int)
    timepoint_df = subset(pretwitch_df, :time => _time -> _time .== time)
    return Dict(timepoint_df.cell .=> dataframe_to_point3f_array(timepoint_df))
end



function pretwitch_over_time(pretwitch_df)
    f = Figure()
    ax = Axis3(f[1, 1], aspect=:data)
    s = Slider(f[2, 1], range=0:360, startvalue=0)
    pt = get_pretwitch_points_at_time(pretwitch_df, 360) |> values |> collect
    sc = Makie.meshscatter!(ax, pt, markersize=5, color=:blue)
    #pt = get_pretwitch_points_at_time(pretwitch_df, 0)
    #sc.positions[] = pt
    on(s.value) do v
        pt = get_pretwitch_points_at_time(pretwitch_df, v) |> values |> collect
        sc.positions[] = pt
    end
    return f
end

function get_pretwitch_annotation_points_at_time(time::Int; pretwitch_df=get_pretwitch_df())
    timepoint_df = subset(pretwitch_df, :time => _time -> _time .== time)
    annotation_points = Dict{String,Point3f}()
    for row in eachrow(timepoint_df)
        cell = strip(row.cell)
        annotation_points[cell] = Point3f(row.x, row.y, row.z)
    end
    return annotation_points
end

function pretwitch_seamcells_over_time(;
    pretwitch_df=get_pretwitch_df(),
    avg_models=nothing,
    average_annotations_dict=nothing,
    left_seam_cells=["a0L"; left_seam_cells],
    right_seam_cells=["a0R"; right_seam_cells],
    seam_cell_midpoints=["a0M"; seam_cell_midpoints],
    voxel_size=0.1625
)
    f = Figure()
    ax = Axis3(f[1, 1], aspect=:data,
        xlabel="L - R (μm)",
        ylabel="A - P (μm)",
        zlabel="D - V (μm)",
        titlealign=:left
    )
    xlims!(ax, -15, 15)
    zlims!(ax, -15, 15)
    #ylims!(ax, -75,125)
    ylims!(ax, 0, 200)
    s = Slider(f[2, 1], range=20:380, startvalue=0)
    seam_cells = vcat(left_seam_cells, right_seam_cells)

    if !isnothing(average_annotations_dict)
        posttwitch_annotations_over_time = map(1:201) do t
            vcat(map(values(average_annotations_dict)) do v
                v.positions[t]
            end...)
        end
    end

    N_prewitch_timepoints = 361
    points_per_seam_cell = map(seam_cells) do seam_cell
        # returns 361 points for the seam cells which are time points 20 to 380 in pretwitch_df
        if seam_cell in ["a0L", "a0R"]
            points = repeat([Point3f(NaN, NaN, NaN)], N_prewitch_timepoints)
        else
            points = get_full_lineage_seam_cell_points(seam_cell; pretwitch_df=pretwitch_df)
        end
        #=points .= map(points) do p
            # Point3(p[3], p[1], p[1])
        end
        =#
        points .= ([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0],) .* points
        N_pretwitch_points = length(points)
        if !isnothing(avg_models)
            # avg_model seam cells are all marked L for some reason
            left_seam_cell = replace(seam_cell, 'R' => 'L')
            seam_cell_index = findfirst(isequal(left_seam_cell), avg_models[1].names)
            seam_cell_index = (seam_cell_index + 1) ÷ 2
            if seam_cell == left_seam_cell
                seam_cell_index += 11
            end
            avg_models_points = map(avg_models) do avg_model
                model_pts = seam_cell_pts(avg_model, 2)
                local p = model_pts[seam_cell_index]
                Point3f(p[1], p[3], p[2])
            end
            append!(points, avg_models_points)
        end
        return seam_cell => points .* voxel_size
    end |> Dict{String,Vector{Point3f}}

    last_pretwitch_H2M = (points_per_seam_cell["H2L"][N_prewitch_timepoints] +
                          points_per_seam_cell["H2R"][N_prewitch_timepoints]) / 2
    first_posttwitch_H2M = last_pretwitch_H2M

    if !isnothing(avg_models)
        @info "Loaded avg models seam cell points"
        # s.range[] = 20:(length(first(values(points_per_seam_cell)))-1)+20
        s.range[] = [20:380; LinRange(380, 381, 11)[2:end-1]; LinRange(381, 751, length(avg_models))]  # pretwitch + avg models
        first_posttwitch_H2M = (points_per_seam_cell["H2L"][N_prewitch_timepoints+1] +
                                points_per_seam_cell["H2R"][N_prewitch_timepoints+1]) / 2
    end

    #pretwitch_translation = first_posttwitch_H2M - last_pretwitch_H2M
    pretwitch_translation = last_pretwitch_H2M - first_posttwitch_H2M

    left_pt = map(left_seam_cells) do seam_cell
        points_per_seam_cell[seam_cell][end]
    end
    right_pt = map(right_seam_cells) do seam_cell
        points_per_seam_cell[seam_cell][end]
    end
    mid_pt = (left_pt .+ right_pt) / 2

    # Translate all points such that H2M is at the origin
    H2_index = findfirst(isequal("H2M"), seam_cell_midpoints)
    H2_mid_pt = mid_pt[H2_index]
    #=
    left_pt .-= H2_mid_pt
    right_pt .-= H2_mid_pt
    mid_pt .-= H2_mid_pt
    =#

    annotation_pts = collect(values(get_pretwitch_annotation_points_at_time(N_prewitch_timepoints; pretwitch_df)))
    annotation_pts .= ([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0],) .* annotation_pts .* voxel_size
    annotation_pts .-= pretwitch_translation

    N = length(left_seam_cells)

    faces = vcat(
        GeometryBasics.NgonFace{3,Int}[(i, i + 1, i + N) for i in 1:N-1],
        GeometryBasics.NgonFace{3,Int}[(i + 1, i + N, i + 1 + N) for i in 1:N-1]
    )
    left_mesh = GeometryBasics.Mesh(
        vcat(left_pt, mid_pt),
        faces
    ) |> Observable
    right_mesh = GeometryBasics.Mesh(
        vcat(right_pt, mid_pt),
        faces
    ) |> Observable

    Makie.mesh!(ax, left_mesh, color=:red, shading=false, alpha=0.5, transparency=true)
    Makie.mesh!(ax, right_mesh, color=:green, shading=false, alpha=0.5, transparency=true)
    left_sc = Makie.meshscatter!(ax, left_pt, markersize=1 * voxel_size, color=:red)
    right_sc = Makie.meshscatter!(ax, right_pt, markersize=1 * voxel_size, color=:green)
    mid_sc = Makie.meshscatter!(ax, mid_pt, markersize=2 * voxel_size, color=:magenta)
    annotation_sc = Makie.meshscatter!(ax, annotation_pts, markersize=2 * voxel_size, color=:cyan)
    left_lines = lines!(ax, left_pt, color=:red)
    right_lines = lines!(ax, right_pt, color=:green)
    mid_lines = lines!(ax, mid_pt, color=:magenta)
    left_text = text!(ax, left_pt; text=left_seam_cells, align=(:center, :bottom))
    right_text = text!(ax, right_pt; text=right_seam_cells, align=(:center, :bottom))
    mid_text = text!(ax, mid_pt; text=seam_cell_midpoints, align=(:center, :bottom))
    on(s.value) do t
        ax.title[] = "Timepoint: $(@sprintf("%06.2f", t)) mpfc $(t < 380 ? "(pretwitch)" : t < 381 ? "(straightening)" : "(posttwitch)")"
        v = if t < 381
            floor(Int, t - 20 + 1)  # pretwitch timepoints start at 20
        else
            round(Int, (t - 381) / ((751 - 381) / 200) + 381) - 20 + 1  # avg models start after pretwitch timepoints
        end
        @info "" v
        # Get points at time v
        left_pt = map(left_seam_cells) do seam_cell
            points_per_seam_cell[seam_cell][v]
        end
        right_pt = map(right_seam_cells) do seam_cell
            points_per_seam_cell[seam_cell][v]
        end

        if t > 380 && t < 381 && !isnothing(avg_models)
            L = t - 380
            @info "Interpolating avg models seam cell points at L=$L"
            function seam_cell_interp_point(seam_cell)
                last_pretwitch_pt = points_per_seam_cell[seam_cell][N_prewitch_timepoints]
                last_pretwitch_pt -= pretwitch_translation
                first_posttwitch_pt = points_per_seam_cell[seam_cell][N_prewitch_timepoints+1]
                last_pretwitch_pt * (1 - L) + first_posttwitch_pt * L
            end
            left_pt = map(seam_cell_interp_point, left_seam_cells)
            right_pt = map(seam_cell_interp_point, right_seam_cells)
        end

        mid_pt = (left_pt .+ right_pt) / 2

        # Translate all points such that H2M is at the origin
        H2_mid_pt = mid_pt[H2_index]

        #  Translate point
        #=
        left_pt .-= H2_mid_pt
        right_pt .-= H2_mid_pt
        mid_pt .-= H2_mid_pt
        =#
        if t < 380
            left_pt .-= pretwitch_translation
            right_pt .-= pretwitch_translation
            mid_pt .-= pretwitch_translation
            println("Applied pretwitch translation: ", pretwitch_translation)
            println("H2_mid_pt: ", H2_mid_pt)
        end

        #left_mesh[].vertices = vcat(left_pt, mid_pt)
        #right_mesh[].vertices = vcat(right_pt, mid_pt)
        left_mesh[] = GeometryBasics.Mesh(
            vcat(left_pt, mid_pt),
            faces
        )
        right_mesh[] = GeometryBasics.Mesh(
            vcat(right_pt, mid_pt),
            faces
        )
        #notify(left_mesh)
        #notify(right_mesh)

        # Update scatter positions
        left_sc.positions[] = left_pt
        right_sc.positions[] = right_pt
        mid_sc.positions[] = mid_pt

        # Update annotation positions
        if isnothing(average_annotations_dict) || t < 380
            annotation_pts = collect(values(get_pretwitch_annotation_points_at_time(
                min(v, N_prewitch_timepoints); pretwitch_df
            )))
            annotation_pts .= ([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0],) .* annotation_pts .* voxel_size
            annotation_pts .-= pretwitch_translation
        elseif !isnothing(average_annotations_dict) && t >= 381
            annotation_pts = posttwitch_annotations_over_time[v-361]
        end
        annotation_sc.positions[] = annotation_pts

        #println("Annotation pts: ", annotation_pts)
        #println(left_pt)

        # Update text positions
        left_text.positions[] = left_pt
        right_text.positions[] = right_pt
        mid_text.positions[] = mid_pt

        # Update line positions
        left_lines.positions[] = left_pt
        right_lines.positions[] = right_pt
        mid_lines.positions[] = mid_pt
    end
    return f
end

function record_pretwitch_seamcells_over_time(;
    pretwitch_df=get_pretwitch_df(),
    avg_models=nothing,
    average_annotations_dict=nothing,
)
    f = pretwitch_seamcells_over_time(; pretwitch_df=pretwitch_df, avg_models=avg_models, average_annotations_dict)
    slider = f.content[2]
    display(f)
    readline()
    record(f, "pretwitch_seamcells_over_time.mp4", 1:length(f.content[2].range[]); framerate=30) do i
        set_close_to!(slider, slider.range[][i])
    end
end

function get_full_lineage_seam_cell_points(seam_cell::String; pretwitch_df=get_pretwitch_df())
    lineage = seam_cell_to_lineage_map[seam_cell]
    df = get_full_lineage_df(lineage; pretwitch_df=pretwitch_df)
    return dataframe_to_point3f_array(df)
end

function get_lattice_points(; pretwitch_df=get_pretwitch_df())
    seam_cells = vcat(left_seam_cells, right_seam_cells)
    points = Dict{String,Vector{Point3f}}()
    for seam_cell in seam_cells
        points[seam_cell] = get_full_lineage_seam_cell_points(seam_cell; pretwitch_df=pretwitch_df)
    end
    return points
end

function get_lattice_points_by_time(; pretwitch_df=get_pretwitch_df())
    seam_cells = vcat(left_seam_cells, right_seam_cells)
    lattice_points = get_lattice_points(; pretwitch_df)
    lattice_points_by_time = Vector{Vector{Point3f}}(undef, 361)
    for t in 1:361
        lattice_points_by_time[t] = map(seam_cells) do seam_cell
            lattice_points[seam_cell][t]
        end
    end
    return lattice_points_by_time
end

function plot_lattice_at_time(
    time::Int,
    origin,
    rotation=I(3);
    pretwitch_df=get_pretwitch_df(),
    lattice_points_by_time=get_lattice_points_by_time(; pretwitch_df=pretwitch_df),
    voxel_size=0.1625
)
    pretwitch_t = get_pretwitch_points_at_time(pretwitch_df, time) |> values |> collect
    pretwitch_t_points = pretwitch_t .- origin
    pretwitch_t_points .= (rotation,) .* pretwitch_t_points
    pretwitch_t_points = Point3f.(pretwitch_t_points)
    pretwitch_t_points .*= voxel_size
    s = scatter(pretwitch_t_points)
    left_lattice_points = lattice_points_by_time[time+1][1:length(left_seam_cells)] .- origin
    left_lattice_points .= (rotation,) .* left_lattice_points
    left_lattice_points = Point3f.(left_lattice_points)
    left_lattice_points .*= voxel_size
    right_lattice_points = lattice_points_by_time[time+1][length(left_seam_cells)+1:end] .- origin
    right_lattice_points .= (rotation,) .* right_lattice_points
    right_lattice_points = Point3f.(right_lattice_points)
    right_lattice_points .*= voxel_size

    mid_lattice_points = (left_lattice_points .+ right_lattice_points) / 2

    lines!(left_lattice_points, color=:red)
    text!(left_lattice_points; text=left_seam_cells)
    lines!(right_lattice_points, color=:green)
    text!(right_lattice_points; text=right_seam_cells)
    lines!(mid_lattice_points, color=:magenta)
    text!(mid_lattice_points; text=seam_cell_midpoints)
    return s
end

# Rotate vector a to align with vector b
function get_rotation_matrix(a, b)
    u_a = normalize(a)
    u_b = normalize(b)
    v = cross(u_a, u_b) # Rotation axis (scaled by sine)
    c = dot(u_a, u_b)   # Cosine of the angle
    s = norm(v)         # Sine of the angle
    if s < 1e-10
        R = c > 0 ? I(3) : -I(3) # Identity if same, -Identity if opposite
    else
        # Skew-symmetric cross-product matrix of v
        vx = [0 -v[3] v[2];
            v[3] 0 -v[1];
            -v[2] v[1] 0]

        # Rodrigues' Rotation Formula
        R = I(3) + vx + (vx^2 * (1 - c) / (s^2))
    end
    return R
end

function get_mid_points_by_time(; pretwitch_df=get_pretwitch_df())
    lattice_points_by_time = get_lattice_points_by_time(; pretwitch_df)
    mid_points_by_time = map(lattice_points_by_time) do lattice_points
        (@view(lattice_points[1:10]) .+ @view(lattice_points[11:20])) / 2
    end
end

# Standard rotation matrices used by interactive plotting helpers below.
# Wrapped as a function so we don't compute (and cache as globals) on
# every package load.
function get_standard_rotations()
    return (;
        R_x_to_y = get_rotation_matrix([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        R_x_to_z = get_rotation_matrix([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        R_rotate_90_about_x = get_rotation_matrix([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
    )
end

# Compute H0 / H2 lineage points and their midpoints from the pretwitch
# coords. Returns a NamedTuple. Lazy — only called when needed; no
# top-level CSV read at module load.
function get_H_points(; pretwitch_df = get_pretwitch_df())
    H0L_points = get_full_lineage_df(seam_cell_to_lineage_map["H0L"]; pretwitch_df) |> dataframe_to_point3f_array
    H0R_points = get_full_lineage_df(seam_cell_to_lineage_map["H0R"]; pretwitch_df) |> dataframe_to_point3f_array
    H2L_points = get_full_lineage_df(seam_cell_to_lineage_map["H2L"]; pretwitch_df) |> dataframe_to_point3f_array
    H2R_points = get_full_lineage_df(seam_cell_to_lineage_map["H2R"]; pretwitch_df) |> dataframe_to_point3f_array
    return (;
        H0L_points, H0R_points,
        H0_midpoints = (H0L_points .+ H0R_points) ./ 2,
        H2L_points, H2R_points,
        H2_midpoints = (H2L_points .+ H2R_points) ./ 2,
    )
end

function get_annotation_name_translation_df(path::AbstractString = _NAMING_CORRELATIONS_CSV)
    return CSV.read(path, DataFrame)
end

function get_color_code_df(path::AbstractString = _COLOR_CODE_CSV)
    return CSV.read(path, DataFrame)
end

function best_annotation_name_match(annotation_name; annotation_name_translation_df=get_annotation_name_translation_df())

end

function audit_annotation_names(flattened_datasets; annotation_name_translation_df=get_annotation_name_translation_df())
    unknown_df = DataFrame(cell_key_path=String[], annotation=String[], best_guess=String[])
    foreach(flattened_datasets) do dataset
        foreach(values(dataset.cell_key.mapping)) do v
            in_ryan = v ∈ annotation_name_translation_df.var"Positional Model Cell Name"
            if !in_ryan
                processed_v = lowercase(strip(v))
                best_guess = ""
                foreach(annotation_name_translation_df.var"Positional Model Cell Name") do ryan_name
                    try
                        processed_ryan_name = lowercase(strip(ryan_name))
                        if processed_v == processed_ryan_name
                            best_guess = ryan_name
                        end
                    catch err
                        @error "Error processing annotation name" v ryan_name err
                    end
                end
                cell_key_path = joinpath(dirname(dataset.path), "CellKey.csv")
                @info "" cell_key_path v best_guess
                push!(unknown_df, (cell_key_path, v, best_guess))
            end
        end
    end
    return unknown_df
end

function get_pretwitch_explicit_df(
    pretwitch_df=get_pretwitch_df();
    use_micrometers=true,
    voxel_size=0.1625,
    avg_models=nothing
)
    last_pretwitch_time_df = subset(pretwitch_df, :time => ByRow(==(360)))
    last_pretwitch_H2L_df = subset(last_pretwitch_time_df,
        :cell => ByRow(==(seam_cell_to_lineage_map["H2L"]))
    )
    last_pretwitch_H2R_df = subset(last_pretwitch_time_df,
        :cell => ByRow(==(seam_cell_to_lineage_map["H2R"]))
    )
    last_pretwitch_H2M = (
        Point3f(last_pretwitch_H2L_df.x[1], last_pretwitch_H2L_df.y[1], last_pretwitch_H2L_df.z[1]) +
        Point3f(last_pretwitch_H2R_df.x[1], last_pretwitch_H2R_df.y[1], last_pretwitch_H2R_df.z[1])
    ) / 2

    if isnothing(avg_models)
        @info "No avg models provided, using pretwitch H2M as reference point"
        first_posttwitch_H2M = Point3f(0.0, 0.0, 0.0)
    else
        @info "Avg models provided, using first posttwitch H2M as reference point"
        local pts = seam_cell_pts(avg_models[1], 2)
        seam_cell_index = findfirst(isequal("H2L"), avg_models[1].names)
        seam_cell_index = (seam_cell_index + 1) ÷ 2
        H2L_index = seam_cell_index + 11
        H2R_index = seam_cell_index
        first_posttwitch_H2M = (pts[H2L_index] + pts[H2R_index]) / 2
    end

    pretwitch_df = copy(pretwitch_df)
    # MPFC
    pretwitch_df.time .+= 20 
    # AP
    pretwitch_df.x .-= last_pretwitch_H2M[1] - first_posttwitch_H2M[3]
    # DV
    pretwitch_df.y .-= last_pretwitch_H2M[2] - first_posttwitch_H2M[2]
    # LR
    pretwitch_df.z .-= last_pretwitch_H2M[3] - first_posttwitch_H2M[1]
    if use_micrometers
        return DataFrame(
            lineage_name = pretwitch_df.cell,
            minutes_post_first_cleavage = pretwitch_df.time,
            LR_micrometers = pretwitch_df.z .* voxel_size,
            DV_micrometers = pretwitch_df.y .* voxel_size,
            AP_micrometers = pretwitch_df.x .* voxel_size
        )
    else
        return DataFrame(
            lineage_name = pretwitch_df.cell,
            minutes_post_first_cleavage = pretwitch_df.time,
            LR_voxels = pretwitch_df.z,
            DV_voxels = pretwitch_df.y,
            AP_voxels = pretwitch_df.x
        )
    end
end