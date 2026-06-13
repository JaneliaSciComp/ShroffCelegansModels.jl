using WGLMakie
using Bonito
using Makie
using Printf
using GeometryBasics
using ShroffCelegansModels.CSV
using ShroffCelegansModels.DataFrames: DataFrame, eachrow

# ─────────────────────────────────────────────────────────────────────────────
# Server-less interactive HTML export of the average-annotation meshscatter.
#
# A static export (Bonito.export_static) has NO Julia backend, so Makie `lift`
# closures never run — a Makie.Slider would leave the scene frozen. Following
# mkitti/BrownianMotionDemo.jl, we instead:
#   • build the scene with STATIC positions (last timepoint),
#   • drive time entirely in the browser with a native HTML `Bonito.Slider`
#     (requires `use_html_widgets = true`),
#   • on slider change, overwrite each meshscatter's GPU instance-position buffer
#     (`positions_transformed_f32c`) via `onjs`, with every timepoint's flat Float32
#     positions embedded at export time,
#   • render the HPF label as a DOM overlay (animates via onjs AND stays visible
#     through zoom, unlike in-scene text!).
#
# Layout follows the demo: a fixed-size Figure inside a centered, max-width flex
# column — NOT resize_to=:body / 100vw, which overflows into a horizontal scrollbar.
# This script is intentionally SELF-CONTAINED (it does not share the movie builder
# meshscatter_average_simple.jl, which uses a Makie.Slider + lift for CairoMakie/
# GLMakie recording).
# ─────────────────────────────────────────────────────────────────────────────

# WGLMakie's internal flat Float32 instance-position GPU buffer for meshscatter.
# Verified present in the pinned WGLMakie (plot-primitives.jl). If a WGLMakie bump
# renames it, this is the single line to update.
const POS_BUFFER_KEY = "positions_transformed_f32c"

function _load_colors_dict()
    colors_csv = joinpath(@__DIR__, "..", "colors.csv")
    colors_dict = Dict{String, RGBf}()
    colors_df = CSV.read(colors_csv, DataFrame)
    for row in eachrow(colors_df)
        colors_dict[lowercase(row.Annotation)] = RGBf(row.R/255, row.G/255, row.B/255)
    end
    return colors_dict
end

const _LEGEND_COLORS = RGBf[
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
const _LEGEND_TEXT = String[
    "Neurons", "Neuroblast", "Muscle", "Intestine", "Nerve Ring",
    "Glial", "Seam Cell", "Other Hypodermal", "Pharyngeal",
]

# avg_dict: Dict from load_latest_average_annotations()
# output_path: destination .html file path
# view: :yz (default), :xz, or :xy
# figure_size / max_width_px: the fixed figure size and the centered card cap.
function export_meshscatter_static(average_annotations_dict;
        output_path::String,
        view::Symbol = :yz,
        figure_size = (1536, 560),
        max_width_px::Int = 1536)
    # resize_to=:parent makes WGLMakie size the canvas to its parent element
    # CLIENT-SIDE at setup (initialize_canvas_size), so it fills the container on
    # first paint instead of coming up small until the first interaction. The parent
    # (fig_box) must have a real height — supplied via aspect-ratio below.
    WGLMakie.activate!(use_html_widgets = true, resize_to = :parent)

    colors_dict = _load_colors_dict()
    get_color(a) = get(colors_dict, lowercase(a), RGBAf(1, 1, 1, 1))

    coordinates = values(average_annotations_dict)
    time_points = axes(first(coordinates).positions, 1)
    T = length(time_points)
    t0 = first(time_points)
    t_last = last(time_points)

    # Initial HPF string (slider starts at the last timepoint). JS mirrors this formula.
    _init_total = (t_last - 1 + 21) / (T - 1) * 370
    _init_hpf = "hpf = $(6 + floor(Int, _init_total / 60)):$(lpad(floor(Int, mod(_init_total, 60)), 2, '0'))"

    label_offset = 8
    scalebar_size_um = 10
    scalebar_y_offset = 180

    app = App(; title = "C. elegans Average Annotations") do session::Session
        with_theme(theme_black()) do
            fig = Figure(size = figure_size)
            ax = LScene(fig[1, 1]; show_axis = false)

            # Native HTML slider (the "scrubber") — drives the animation client-side.
            sg = Bonito.Slider(time_points, value = t_last)

            # Static meshscatter per group; collect plot handle + per-timepoint flat
            # Float32 positions (interleaved x,y,z) for the onjs buffer update.
            # NOTE: explicit flatten — reinterpret(Float32, ::Vector{Point3{Float64}})
            # is invalid (Float64→Float32 width mismatch). Rendered unfiltered so the
            # point count is constant across timepoints (required: fixed-size GPU buffer).
            plots = Any[]
            group_positions = Vector{Vector{Float32}}[]
            # uuid → per-point annotation names, for the client-side hover tooltip.
            # js_uuid(plot) == string(objectid(plot)) matches the JS `plot_uuid`
            # returned by WGL.pick_closest.
            names_by_uuid = Dict{String, Vector{String}}()
            for v in coordinates
                cols = get_color.(v.annotations)
                p = meshscatter!(ax, v.positions[t_last]; markersize = 1.0, color = cols,
                    inspectable = true)
                push!(plots, p)
                names_by_uuid[WGLMakie.js_uuid(p)] = String.(v.annotations)
                push!(group_positions,
                      [Float32[c for pt in v.positions[t] for c in pt] for t in time_points])
            end

            # Animate positions client-side: on slider change, overwrite each group's
            # GPU position buffer with that timepoint's embedded flat Float32 array.
            for (p, perT) in zip(plots, group_positions)
                onjs(session, sg.value, js"""(val) => {
                    const data = $(perT);
                    const t = Math.round(val) - $(t0);    // slider value -> 0-based index
                    if (t < 0 || t >= data.length) return;
                    $(p).then(plots => {
                        plots[0].plot_object.update([[$(POS_BUFFER_KEY), new Float32Array(data[t])]]);
                    });
                }""")
            end

            # Scalebar LINE stays in data space (overdraw=true) so its length tracks
            # zoom correctly. The LABEL is a DOM overlay (below) — in-scene text! in
            # WGLMakie 3D is orientation-dependent (only readable from some angles),
            # whereas a DOM label is always visible.
            scale_str, scalebar = if view == :xy
                ("$(scalebar_size_um ÷ 2) μm",
                 Point3f[
                     [label_offset+1,                        0, -label_offset-1],
                     [label_offset+1 - scalebar_size_um/2,   0, -label_offset-1],
                 ])
            else
                ("$scalebar_size_um μm",
                 Point3f[
                     [label_offset+1, scalebar_y_offset,                    -label_offset-1],
                     [label_offset+1, scalebar_y_offset + scalebar_size_um, -label_offset-1],
                 ])
            end
            lines!(ax, scalebar, color = :white, linewidth = 5, overdraw = true)

            # Legend: a DOM overlay (built below) rather than a Makie Legend. A Makie
            # legend is baked into the WGLMakie scene at export time and can't reflow on
            # window resize in a static file; an HTML flex-wrap row reflows natively.

            cc = Camera3D(ax.scene;
                projectiontype = Makie.Orthographic,
                lookat = Vec3d(0, 90, 0), eyeposition = Vec3d(60, 90, 0))
            zoom!(ax.scene, 0.30)
            if view == :yz
                update_cam!(ax.scene, cc, 0, 0)
            elseif view == :xz
                update_cam!(ax.scene, cc, 0, π/2)
            elseif view == :xy
                update_cam!(ax.scene, cc, π/2, 0)
            else
                error("view must be :yz, :xz, or :xy; got $view")
            end

            # HPF label as an absolutely-positioned DOM overlay over the figure box:
            # animates via onjs and never disappears under the 3D camera.
            hpf = DOM.div(_init_hpf; id = "hpf-label", style = Styles(CSS(
                "position" => "absolute", "top" => "10px", "left" => "14px",
                "color" => "white", "font-size" => "20px",
                "font-family" => "sans-serif", "z-index" => "10",
                "pointer-events" => "none")))

            # Scalebar label as a DOM overlay (always visible; pairs with the in-scene
            # line). Top-right so it clears the bottom DOM legend.
            scale_label = DOM.div("Scale Bar: $scale_str"; id = "scale-label", style = Styles(CSS(
                "position" => "absolute", "top" => "10px", "right" => "14px",
                "color" => "white", "font-size" => "18px",
                "font-family" => "sans-serif", "z-index" => "10",
                "pointer-events" => "none")))

            # Legend as a DOM overlay across the bottom: an HTML flex-wrap row reflows
            # natively as the window/figure width changes (a Makie legend can't, since
            # it's frozen into the scene at export time). pointer-events:none so it never
            # blocks camera interaction.
            legend_items = map(zip(_LEGEND_COLORS, _LEGEND_TEXT)) do (c, label)
                rgb = "rgb($(round(Int, 255*c.r)),$(round(Int, 255*c.g)),$(round(Int, 255*c.b)))"
                DOM.span(
                    DOM.span(""; style = Styles(CSS(
                        "display" => "inline-block", "width" => "13px", "height" => "13px",
                        "border-radius" => "50%", "background-color" => rgb,
                        "margin-right" => "6px", "flex" => "0 0 auto"))),
                    DOM.span(label);
                    style = Styles(CSS("display" => "inline-flex", "align-items" => "center",
                        "white-space" => "nowrap")))
            end
            legend = DOM.div(legend_items...; id = "legend", style = Styles(CSS(
                "position" => "absolute", "bottom" => "6px", "left" => "0", "right" => "0",
                "display" => "flex", "flex-wrap" => "wrap", "justify-content" => "center",
                "align-items" => "center", "gap" => "4px 16px",
                "padding" => "0 12px", "box-sizing" => "border-box",
                "color" => "white", "font-family" => "sans-serif", "font-size" => "15px",
                "z-index" => "10", "pointer-events" => "none")))
            onjs(session, sg.value, js"""(val) => {
                const T = $(T);
                const t = Math.round(val);
                const total = (t - 1 + 21) / (T - 1) * 370;
                const hours = 6 + Math.floor(total / 60);
                const minutes = Math.floor(total % 60);
                const el = document.getElementById("hpf-label");
                if (el) el.textContent = "hpf = " + hours + ":" + String(minutes).padStart(2, "0");
            }""")

            # DataInspector-style hover tooltip. The standard Makie DataInspector
            # computes its label via a Julia callback on mouse events, which can't run
            # in a server-less export. Instead we do the picking + label entirely
            # client-side: on mousemove, WGL.pick_closest returns [plot_uuid, index];
            # we look the annotation name up in the embedded uuid→names map and show a
            # DOM tooltip at the cursor. position:fixed + clientX/Y keeps placement
            # correct regardless of ancestor positioning/scroll.
            tooltip = DOM.div(""; id = "inspector-tooltip", style = Styles(CSS(
                "position" => "fixed", "display" => "none", "z-index" => "20",
                "background" => "rgba(0,0,0,0.82)", "color" => "white",
                "padding" => "4px 8px", "border-radius" => "4px",
                "font-family" => "sans-serif", "font-size" => "14px",
                "pointer-events" => "none", "white-space" => "nowrap")))
            Bonito.evaljs(session, js"""
                Promise.all([$(WGLMakie.WGL), $(ax.scene)]).then(([WGL, scene]) => {
                    if (!scene || !scene.screen) { return; }
                    const canvas = scene.screen.canvas;
                    const lookup = $(names_by_uuid);
                    const tip = $(tooltip);
                    const POS_KEY = $(POS_BUFFER_KEY);
                    canvas.addEventListener("mousemove", (event) => {
                        const xy = WGL.events2unitless(scene.screen, event);
                        // 1x1 pick: only fires when the cursor is actually over a point
                        // (empty area -> picks.length == 0 -> hide), so no nearest-point
                        // tooltip in blank space.
                        const picked = WGL.pick_native(scene, xy[0], xy[1], 1, 1);
                        if (picked) {
                            const picks = picked[1];
                            if (picks.length === 1) {
                                const [plot, index] = picks[0];
                                const names = lookup[plot.plot_uuid];
                                if (names) {
                                    const name = (names[index] !== undefined) ? names[index] : "?";
                                    // Current coordinate straight from the GPU instance
                                    // buffer (reflects the animated timepoint).
                                    let coord = "";
                                    const attr = plot.geometry.attributes[POS_KEY];
                                    if (attr && attr.array) {
                                        const a = attr.array;
                                        const x = a[index*3], y = a[index*3+1], z = a[index*3+2];
                                        coord = " (" + x.toFixed(1) + ", " + y.toFixed(1) + ", " + z.toFixed(1) + ")";
                                    }
                                    tip.innerText = name + coord;
                                    tip.style.left = (event.clientX + 12) + "px";
                                    tip.style.top  = (event.clientY + 12) + "px";
                                    tip.style.display = "block";
                                    return;
                                }
                            }
                        }
                        tip.style.display = "none";
                    });
                    canvas.addEventListener("mouseleave", () => { tip.style.display = "none"; });

                    // Camera reconcile kick. In attach_3d_camera, the makie projection
                    // is only rebuilt (update_matrices) from the OrbitControls "change"
                    // handler — which never fires at load, so the plot renders with the
                    // serialized projection until the first scroll/drag. The correct fix
                    // is OrbitControls.update(): the SAME call scroll/drag make. It
                    // reconciles the controls' spherical state with the camera and THEN
                    // dispatches "change" (a bare dispatchEvent skips the reconciliation
                    // and mis-frames the view). allow_update() is true in a server-less
                    // export (no Julia), so this runs client-side. Fire across a few
                    // frames to catch the settled canvas resolution.
                    const recompute = () => {
                        if (scene.orbitcontrols) { scene.orbitcontrols.update(); }
                    };
                    requestAnimationFrame(recompute);
                    setTimeout(recompute, 50);
                    setTimeout(recompute, 250);
                    setTimeout(recompute, 600);
                });
            """)

            # Centered, max-width flex column (demo pattern): no element exceeds the
            # viewport width, so there is no horizontal scrollbar on any screen size.
            # Anchor the overlays to the figure box (position:relative) so they sit in
            # the figure's corners rather than the whole card.
            fig_box = DOM.div(hpf, scale_label, legend, tooltip, fig; style = Styles(CSS(
                "position" => "relative", "width" => "100%",
                "aspect-ratio" => "$(figure_size[1]) / $(figure_size[2])")))

            # Make the native range input stretch to fill its (flex) container.
            slider_css = DOM.style("input[type=range]{width:100%; accent-color:#38bdf8; height:6px;}")
            card = DOM.div(
                slider_css, fig_box,
                DOM.div(
                    DOM.span("Time"; style = Styles(CSS(
                        "color" => "#cbd5e1", "font-weight" => "700",
                        "font-family" => "sans-serif", "font-size" => "0.95rem"))),
                    DOM.div(sg; style = Styles(CSS("flex" => "1", "width" => "100%")));
                    style = Styles(CSS("width" => "100%", "display" => "flex",
                        "align-items" => "center", "gap" => "0.75rem"))),
                style = Styles(CSS(
                    "position" => "relative", "width" => "100%",
                    "max-width" => "$(max_width_px)px",
                    "display" => "flex", "flex-direction" => "column",
                    "align-items" => "center", "gap" => "0.5rem")))

            # Zero the default <html>/<body> margin and paint them black, otherwise the
            # browser's default body margin shows as a white border around the card.
            body_css = DOM.style("html,body{margin:0;padding:0;background-color:#000;}")

            return DOM.div(body_css, card; style = Styles(CSS(
                "min-height" => "100vh", "background-color" => "black",
                "margin" => "0", "padding" => "1rem", "box-sizing" => "border-box",
                "display" => "flex", "flex-direction" => "column",
                "align-items" => "center")))
        end
    end
    Bonito.export_static(output_path, app)
    bytes = isfile(output_path) ? filesize(output_path) : 0
    @info "Static export written" path=output_path size_mb=round(bytes/1024^2; digits=1)
end
