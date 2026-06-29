"""
    CommonScenes

Shared precompile coverage for the apps whose own code can't run at build time
(the launch-script apps — see [`ShowAverageAnnotations`](@ref)). Those apps
eagerly read PVC data at load time, so they can't be exercised directly. But
their cold start is dominated by the same Makie → WGLMakie → Bonito
render/serialize compilation, which is shared process-wide and stored in this
package's image. This submodule's `@compile_workload` renders one representative
3-D scene exercising the *union* of plot primitives + widgets those apps use
(`mesh!`, `meshscatter!`, `lines!`, `scatter!`, `text!`, `vlines!`, `LScene`,
`Axis3`, `Axis`, `Toggle`, `Menu`, `Slider`, `Button`, `Label`, `DataInspector`)
and serializes it via `Bonito.export_static`. Every container that does
`using ShroffCelegansModelsWebInterface` loads this cache.

Gated by the same `precompile_workload` Preference as the per-app workloads.
"""
module CommonScenes

using WGLMakie
using Bonito
using GeometryBasics
using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference

"Build a Figure touching the union of primitives the launch-script apps render."
function _representative_figure()
    fig = Figure(size = (1200, 800))

    # 3-D LScene: translucent colored mesh + meshscatter + lines + scatter + text.
    ls = LScene(fig[1, 1]; show_axis = false)
    rectmesh = GeometryBasics.mesh(Rect3f(Vec3f(-1), Vec3f(2)))
    nverts = length(GeometryBasics.coordinates(rectmesh))
    mesh!(ls, rectmesh; color = collect(range(0, 1; length = nverts)),
        colorrange = (0, 1), transparency = true, alpha = 0.3, inspectable = false)
    pts = [Point3f(cos(t), sin(t), 0.1t) for t in range(0, 2π; length = 12)]
    meshscatter!(ls, pts; markersize = 0.2, color = :gray, alpha = 1)
    scatter!(ls, pts; color = :red, markersize = 6)
    lines!(ls, pts; color = :white, linewidth = 2)
    text!(ls, Point3f(0, 0, 2); text = "hpf = 12:31", fontsize = 20)

    # 3-D Axis3: another mesh + scatter (Axis3 codegen differs from LScene).
    ax3 = Axis3(fig[1, 2])
    mesh!(ax3, rectmesh; color = :blue, transparency = true, alpha = 0.2, inspectable = false)
    meshscatter!(ax3, pts; markersize = 0.2, color = :gold)

    # 2-D Axis time series (volumes/lengths over time) + vline.
    ax2 = Axis(fig[2, 1:2])
    xs = collect(1:20)
    lines!(ax2, xs, sin.(xs ./ 3); color = :green)
    vlines!(ax2, [5, 10]; color = :gray)

    # Widgets row.
    grid = GridLayout(fig[3, 1:2])
    grid[1, 1] = Label(fig, "Marker Size")
    grid[1, 2] = Makie.Slider(fig; range = 0.5:0.1:4, startvalue = 1.0)
    grid[1, 3] = Button(fig; label = "XY")
    grid[1, 4] = Toggle(fig; active = true)
    grid[1, 5] = Menu(fig; options = string.(1:5), default = "1")

    DataInspector(fig; backgroundcolor = :black)
    return fig
end

if @load_preference("precompile_workload", true)
@setup_workload begin
    @compile_workload begin
        try
            WGLMakie.activate!()
            app = App(() -> _representative_figure())
            mktempdir() do dir
                export_static(joinpath(dir, "common_scene.html"), app)
            end
        catch err
            @debug "CommonScenes precompile workload skipped" exception = (err, catch_backtrace())
        end
    end
end
end  # precompile_workload preference guard

end # module CommonScenes
