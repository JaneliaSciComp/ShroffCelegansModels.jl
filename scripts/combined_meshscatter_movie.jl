using CairoMakie
using GeometryBasics: Point3f

include(joinpath(@__DIR__, "generate_pipeline_movie.jl"))

# Combined pre-twitch + post-twitch meshscatter movie.
#
# The post-twitch averaged annotations (`avg_dict`, loaded via
# `load_average_annotations`, already in the flipped RL/VD display frame and
# in micrometers) cover mpfc 381→751. This builds a single movie that
# *prepends* the pre-twitch annotation timepoints (mpfc 20→380) read from the
# ryan_data pre-twitch CSV, so one movie spans the whole 20→751 mpfc range in
# the same meshscatter style as the pipeline's :yz / :xz movies.
#
# Cells are shown only while present (per the data): pre-twitch cells are NaN
# during the post-twitch frames and vice-versa. Pre-twitch points are placed
# in the same movie frame as the post-twitch annotations and aligned by
# matching the H2 seam-cell midpoint (H2M) of the last pre-twitch frame to the
# H2M of the first post-twitch frame (the `seam_cells` group in `avg_dict`).

# Raw pre-twitch coords are (x=AP, y=DV, z=LR) in voxels. The movie/display
# frame is (x=LR, y=AP, z=DV) with LR and DV negated (the RL/VD flip), in
# micrometers: (x,y,z) → (-z, x, -y) * voxel_size.
_pretwitch_to_movie(p, voxel_size) = Point3f(-p[3], p[1], -p[2]) * voxel_size

"""
    build_combined_annotations_dict(avg_dict; pretwitch_df, voxel_size)

Return `(combined_dict, frame_mpfc)` where `combined_dict` has the same shape
as a loaded averaged-annotations dict (group → (annotations, positions)) with
the pre-twitch frames prepended, and `frame_mpfc[t]` is the minutes-post-first-
cleavage for frame `t` (used for the hpf label).
"""
function build_combined_annotations_dict(avg_dict;
        pretwitch_df = ShroffCelegansModels.get_pretwitch_df(),
        voxel_size = 0.1625)
    NaN3 = Point3f(NaN, NaN, NaN)
    to_movie(p) = _pretwitch_to_movie(p, voxel_size)

    times = sort(unique(pretwitch_df.time))
    N_pre = length(times)
    N_post = length(first(values(avg_dict)).positions)

    # Per-frame pre-twitch points (raw frame), keyed by lineage cell name.
    pre_points_per_time = [
        ShroffCelegansModels.get_pretwitch_annotation_points_at_time(t; pretwitch_df)
        for t in times
    ]

    # --- Alignment: pre-twitch H2M (last frame) → post-twitch H2M (frame 1). ---
    haskey(avg_dict, "seam_cells") ||
        error("avg_dict has no `seam_cells` group; cannot align pre-twitch to post-twitch")
    seam = avg_dict["seam_cells"]
    h2r_idx = findfirst(==("H2R"), seam.annotations)
    h2l_idx = findfirst(==("H2L"), seam.annotations)
    (isnothing(h2r_idx) || isnothing(h2l_idx)) &&
        error("seam_cells group missing H2R/H2L; got $(seam.annotations)")
    H2M_post = Point3f((seam.positions[1][h2r_idx] .+ seam.positions[1][h2l_idx]) ./ 2)

    h2l_lin = ShroffCelegansModels.seam_cell_to_lineage_map["H2L"]
    h2r_lin = ShroffCelegansModels.seam_cell_to_lineage_map["H2R"]
    last_pre = pre_points_per_time[end]
    (haskey(last_pre, h2l_lin) && haskey(last_pre, h2r_lin)) ||
        error("last pre-twitch frame (t=$(times[end])) missing H2 lineage cells $h2l_lin/$h2r_lin")
    H2M_pre = (to_movie(last_pre[h2l_lin]) + to_movie(last_pre[h2r_lin])) / 2
    translation = H2M_post - H2M_pre

    # --- Pre-twitch group: union of cells over time, NaN where absent. ---
    pre_cells = sort(collect(reduce(union, (Set(keys(p)) for p in pre_points_per_time))))
    pre_positions = Vector{Vector{Point3f}}(undef, N_pre + N_post)
    for (i, pts) in enumerate(pre_points_per_time)
        pre_positions[i] = [haskey(pts, c) ? to_movie(pts[c]) + translation : NaN3
                            for c in pre_cells]
    end
    for i in (N_pre + 1):(N_pre + N_post)
        pre_positions[i] = fill(NaN3, length(pre_cells))
    end

    combined = Dict{String, @NamedTuple{annotations::Vector{String}, positions::Vector{Vector{Point3f}}}}()
    combined["pretwitch"] = (; annotations = pre_cells, positions = pre_positions)

    # --- Post-twitch groups: NaN for the pre-twitch frames, then real data. ---
    for (k, v) in avg_dict
        n = length(v.annotations)
        new_positions = Vector{Vector{Point3f}}(undef, N_pre + N_post)
        for i in 1:N_pre
            new_positions[i] = fill(NaN3, n)
        end
        for i in 1:N_post
            new_positions[N_pre + i] = Point3f.(v.positions[i])
        end
        combined[k] = (; annotations = v.annotations, positions = new_positions)
    end

    # --- Per-frame mpfc for the hpf label. ---
    pre_mpfc = Float64.(times) .+ 20.0                       # raw 0..360 → mpfc 20..380
    post_mpfc = collect(range(381.0, 751.0, length = N_post))
    frame_mpfc = vcat(pre_mpfc, post_mpfc)

    return combined, frame_mpfc
end

"""
    generate_combined_meshscatter_movie(avg_dict; output_path, framerate=48, view=:yz, pretwitch_df)

Render the combined pre+post-twitch meshscatter movie to `output_path`.
Mirrors `generate_meshscatter_movie` but spans both phases.
"""
function generate_combined_meshscatter_movie(avg_dict;
        output_path::String, framerate::Int = 48, view::Symbol = :yz,
        pretwitch_df = ShroffCelegansModels.get_pretwitch_df())
    CairoMakie.activate!()
    combined, frame_mpfc = build_combined_annotations_dict(avg_dict; pretwitch_df)
    raw_path = replace(output_path, ".mp4" => "_raw.mp4")
    with_theme(theme_black()) do
        fig, time_slider, ax = meshscatter_average_simple(combined;
            xy_bounding_radius = -1,
            figure_size = (1920, 360),
            view,
            frame_mpfc)
        time_points = time_slider.range[]
        time_slider.blockscene.visible[] = false
        nframes = length(time_points)
        frame = Ref(0)
        try
            vs = Record(fig, time_points; framerate, px_per_unit = 1, preset = "fast", update = false) do t
                set_close_to!(time_slider, t)
                n = (frame[] += 1)
                @info "combined movie frame  view=$view  frame=$n/$nframes  ($(round(Int, 100n/nframes))%)"
            end
            @info "combined movie frames complete  view=$view  frames=$(frame[])/$nframes  → encoding+crop"
            save(raw_path, vs)
            crop_video(raw_path, output_path; framerate)
        finally
            isfile(raw_path) && rm(raw_path; force = true)
        end
    end
    bytes = isfile(output_path) ? filesize(output_path) : 0
    @info "Combined movie written" path=output_path size_mb=round(bytes/1024^2; digits=1) framerate
end
