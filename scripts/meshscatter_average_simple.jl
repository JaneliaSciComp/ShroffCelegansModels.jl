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
function meshscatter_average_simple(average_annotations_dict; xy_bounding_radius = -1)
    fig = Figure(size = (1920, 1080))
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
    time_slider = Makie.Slider(fig[2, 1], range = time_points, startvalue = last(time_points))

    # HPF label
    time_text = lift(time_slider.value) do t
        total_minutes = (t - 1 + 21) / (length(time_points) - 1) * 370
        hours = 6 + round(Int, total_minutes / 60, RoundDown)
        minutes = round(Int, mod(total_minutes, 60), RoundDown)
        "hpf = $hours:$(@sprintf("%02d", minutes))"
    end
    text!(ax, -label_offset, 0, label_offset; text = time_text, fontsize = _fontsize)

    # Scalebar
    scalebar_size_um = 10
    scalebar_y_offset = 180
    scalebar_text = Observable("$scalebar_size_um μm")
    scalebar = Observable(Point3f[
        [label_offset+1, scalebar_y_offset,                    -label_offset-1],
        [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1]
    ])
    scalebar_label_position = lift(scalebar) do pos
        first(pos)
    end
    text!(ax, scalebar_label_position; text = scalebar_text, fontsize = _fontsize, align = (:left, :bottom))
    lines!(ax, scalebar, color = :white, linewidth = 5)

    alpha_obs = Observable(1.0)

    # Meshscatter for each cell group — lift positions from slider so Bonito
    # can serialize the whole Observable graph for static export.
    for v in coordinates
        colors = get_color.(v.annotations)
        current_pos = if xy_bounding_radius > 0
            lift(time_slider.value) do t
                map(v.positions[t]) do pt
                    pt[1]^2 + pt[3]^2 > xy_bounding_radius^2 ? Point3{Float64}(NaN, NaN, NaN) : pt
                end
            end
        else
            lift(time_slider.value) do t
                v.positions[t]
            end
        end
        meshscatter!(ax, current_pos; markersize = _markersize, color = colors, alpha = alpha_obs)
    end

    # Cell-type color legend
    scatter_legend_pts = Point3f[
        [0,  0, 0],
        [0,  1, 0],
        [0,  2, 0],
        [0,  3, 0],
        [0,  4, 0],
        [0,  5, 0],
        [0,  6, 0],
        [0,  7, 0],
        [0,  8, 0],
    ]
    scatter_legend_pts .*= Point3f(0, 20, 0)
    scatter_legend_pts .+= Point3f(label_offset+3, 0, -label_offset-3)
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
    meshscatter!(ax, scatter_legend_pts; markersize = _markersize, color = legend_colors, alpha = alpha_obs)
    text!(ax, scatter_legend_pts; text = legend_text, align = (:left, :center), offset = (10, 0))

    Camera3D(ax.scene;
        projectiontype = Makie.Orthographic,
        lookat = Vec3d(0, 90, 0),
        eyeposition = Vec3d(60, 90, 0)
    )
    zoom!(ax.scene, 4)

    return fig, time_slider
end
