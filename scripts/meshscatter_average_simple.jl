using Makie
using Makie: throttle
using Printf
using GeometryBasics
using ShroffCelegansModels.CSV
using ShroffCelegansModels.DataFrames: DataFrame, eachrow

function load_colors_dict_simple()
    colors_csv = joinpath(@__DIR__, "..", "colors.csv")
    colors_dict = Dict{String, RGBf}()
    colors_df = CSV.read(colors_csv, DataFrame)
    for row in eachrow(colors_df)
        colors_dict[lowercase(row.Annotation)] = RGBf(row.R/255, row.G/255, row.B/255)
    end
    return colors_dict
end

# Stripped-down meshscatter visualization: time slider only, no record/cosmetic
# controls. Uses lift() so Observable updates are serializable to JavaScript for
# Bonito.export_static and also work with CairoMakie headless recording.
# Returns (fig, time_slider).
function meshscatter_average_simple(average_annotations_dict;
        xy_bounding_radius = -1,
        figure_size = (960, 300),
        show_legend = true,
        view = :yz)
    fig = Figure(size = figure_size)
    ax = LScene(fig[1, 1]; show_axis = false)
    coordinates = values(average_annotations_dict)
    _markersize = Observable(1.0)

    colors_dict = load_colors_dict_simple()
    function get_color(annotation)
        get(colors_dict, lowercase(annotation), RGBAf(1,1,1,1))
    end

    _fontsize = Observable(20)

    label_offset = 8
    time_points = axes(first(coordinates).positions, 1)

    # Create the Slider first so all lift() calls below reference its Observable.
    # Layout: LScene → fig[1,1], Legend → overlaid on fig[1,1], Slider → fig[2,1].
    time_slider = Makie.Slider(fig[2, 1], range = time_points, startvalue = last(time_points))

    # HPF label
    time_text = lift(time_slider.value) do t
        total_minutes = (t - 1 + 21) / (length(time_points) - 1) * 370
        hours = 6 + round(Int, total_minutes / 60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        "hpf = $hours:$(@sprintf("%02d", minutes))"
    end
    text!(ax, -label_offset, 0, label_offset; text = time_text, fontsize = _fontsize)

    # Scalebar — position and label depend on view
    scalebar_size_um = 10
    scalebar_y_offset = 180
    scalebar_text = if view == :xy
        Observable("$(scalebar_size_um÷2) μm")
    else
        Observable("$scalebar_size_um μm")
    end
    scalebar = if view == :xy
        Observable(Point3f[
            [label_offset+1,                      0, -label_offset-1],
            [label_offset+1 - scalebar_size_um/2, 0, -label_offset-1],
        ])
    else
        Observable(Point3f[
            [label_offset+1, scalebar_y_offset,                    -label_offset-1],
            [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1],
        ])
    end
    scalebar_label_position = lift(scalebar) do pos
        first(pos)
    end
    text!(ax, scalebar_label_position; text = scalebar_text, fontsize = _fontsize, align = (:left, :bottom))
    lines!(ax, scalebar, color = :white, linewidth = 5)

    alpha_obs = Observable(1.0)

    # Meshscatter for each cell group.  Uses lift() so the Observable graph is
    # serializable to JavaScript for Bonito.export_static.
    #
    # When xy_bounding_radius > 0, cells outside the cylinder are filtered from
    # BOTH the position and color arrays before passing to meshscatter, so
    # CairoMakie/GLMakie never renders them (no off-screen sentinel overhead).
    # The color array is captured per-group at construction time; the Observable
    # wraps only the filtered (position, color) pair.
    for v in coordinates
        all_colors = get_color.(v.annotations)
        if xy_bounding_radius > 0
            pos_col = lift(time_slider.value) do t
                pts = v.positions[t]
                mask = [pt[1]^2 + pt[3]^2 <= xy_bounding_radius^2 for pt in pts]
                pts[mask], all_colors[mask]
            end
            meshscatter!(ax, lift(pc -> pc[1], pos_col);
                markersize = _markersize,
                color      = lift(pc -> pc[2], pos_col),
                alpha      = alpha_obs)
        else
            current_pos = lift(time_slider.value) do t; v.positions[t]; end
            meshscatter!(ax, current_pos;
                markersize = _markersize, color = all_colors, alpha = alpha_obs)
        end
    end

    legend_colors = RGBf[
        RGBf(185,   0, 255),
        RGBf(  0, 255, 255),
        RGBf(255,   0,   0),
        RGBf(255,  99,  93),
        RGBf(  0,   0, 128),
        RGBf(255, 251,  93),
        RGBf( 16, 185,   0),
        RGBf(  0, 128,   0),
        RGBf( 93, 255, 109),
    ] ./ 255
    legend_text = String[
        "Neurons",
        "Neuroblast",
        "Muscle",
        "Intestine",
        "Nerve Ring",
        "Glial",
        "Seam Cell",
        "Other Hypodermal",
        "Pharyngeal",
    ]
    if show_legend
        # Overlaid on the LScene cell, horizontally aligned with the scalebar (bottom of scene).
        # halign=:left keeps it clear of the scalebar at bottom-right.
        legend_elements = [MarkerElement(color = c, marker = :circle, markersize = 20)
                           for c in legend_colors]
        Legend(fig[1, 1], legend_elements, legend_text;
            orientation     = :horizontal,
            framevisible    = false,
            labelcolor      = :white,
            backgroundcolor = :black,
            labelsize       = 20,
            rowgap          = 0,
            padding         = (4, 4, 4, 4),
            halign          = :center,
            valign          = :bottom,
            tellheight      = false,
            tellwidth       = false,
        )
    end

    cc = Camera3D(ax.scene;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(60, 90, 0)
    )
    zoom!(ax.scene, 0.30)
    # Apply requested view: phi/theta follow update_cam!(scene, cam, phi, theta)
    # YZ side view: eye along +X  → phi=0, theta=0
    # XZ top view:  eye along +Z  → phi=0, theta=π/2
    # XY front view: eye along +Y → phi=π/2, theta=0
    if view == :yz
        update_cam!(ax.scene, cc, 0, 0)
    elseif view == :xz
        update_cam!(ax.scene, cc, 0, π/2)
    elseif view == :xy
        update_cam!(ax.scene, cc, π/2, 0)
    else
        error("view must be :yz, :xz, or :xy; got $view")
    end

    return fig, time_slider, ax
end
