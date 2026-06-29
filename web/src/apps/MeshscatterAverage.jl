"""
    MeshscatterAverage

Submodule for the `meshscatter_average_edited` 3-D meshscatter viewer
(deployment container `meshscatter-average-edited`, port 8580).

The render function `meshscatter_average` (and `load_colors_dict`) are defined
in the shared script `scripts/meshscatter_average_dev.jl`; we `include` it here
so those methods live in this submodule's namespace and become part of this
package's precompile image. The precompile workload at the bottom builds a
figure from tiny synthetic data and serializes it via `Bonito.export_static`,
caching the costly Makie → WGLMakie → Bonito compilation at build time.
"""
module MeshscatterAverage

using WGLMakie
using Bonito
using ShroffCelegansModels
using GeometryBasics: Point3
using ShroffCelegansModels: load_average_annotations, load_latest_average_annotations, load_annotation_cache
using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference

# Defines `meshscatter_average(...)` and `load_colors_dict()` in this module.
# The script lives outside the package source tree; `include` records it as a
# precompile dependency, so edits to it correctly invalidate this package.
include(joinpath(pkgdir(ShroffCelegansModels), "scripts", "meshscatter_average_dev.jl"))

const PROXY_PATH = "meshscatter_average_edited"
const PORT = 8580
const DEFAULT_FILENAME =
    "edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells_2025_02_13.h5"

black_body(fig) = DOM.body(fig, style = Styles(CSS("background-color" => "black")))

"Build the Bonito `App` pair (main + `/nerve_ring`) for a given annotation dict."
function build_apps(average_annotation_dict)
    app = App(; title = "Shroff Lab: C. elegans meshscatter_average") do session::Session
        with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict; session, xy_bounding_radius = sqrt(52)))
        end
    end
    nerve_ring_app = App(; title = "Shroff Lab: C. elegans meshscatter_average/nerve_ring") do
        with_theme(theme_black()) do
            black_body(meshscatter_average(average_annotation_dict; nerve_ring = true))
        end
    end
    return app, nerve_ring_app
end

"Construct (but do not block on) the Bonito `Server` for this app."
function webapp(average_annotation_dict)
    app, nerve_ring_app = build_apps(average_annotation_dict)
    server = Server(app, "0.0.0.0", PORT;
        proxy_url = "https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/$PROXY_PATH/")
    route!(server, "/nerve_ring" => nerve_ring_app)
    return server
end

"Container entry point (called from `web/scripts/web_meshscatter_average_edited.jl`)."
function main()
    @info "Loading annotation cache"
    load_annotation_cache()
    @info "Activating WGLMakie"
    WGLMakie.activate!(; resize_to = :body)
    @info "Launching server!"
    average_annotation_dict = load_latest_average_annotations(default_filename = DEFAULT_FILENAME)
    server = webapp(average_annotation_dict)
    if isinteractive()
        println("Press enter to quit")
        readline()
    else
        wait(server)
    end
end

# --- precompile workload -------------------------------------------------------

"Tiny in-memory stand-in for `load_average_annotations`'s return value."
function _synthetic_annotation_dict(; ntime = 6)
    d = Dict{String, @NamedTuple{annotations::Vector{String}, positions::Vector{Vector{Point3{Float64}}}}}()
    mk(names) = (; annotations = names,
        positions = [[Point3{Float64}(0.1i, 90.0 + j, 0.1j) for j in 1:length(names)] for i in 1:ntime])
    d["DCR6485_RPM1_NU"] = mk(["a1", "a2", "a3", "a4"])  # key the nerve_ring path requires
    d["group_b"] = mk(["b1", "b2", "b3"])
    return d
end

# Toggle the (costly) render workload via a Preference, which — unlike an env
# var — participates in the precompile cache hash, so flipping it reliably
# invalidates the cache. Disable for fast dev iteration with:
#   julia --project=web -e 'using Preferences;
#     set_preferences!("ShroffCelegansModelsWebInterface", "precompile_workload"=>false)'
if @load_preference("precompile_workload", true)
@setup_workload begin
    dict = _synthetic_annotation_dict()
    @compile_workload begin
        # `load_colors_dict()` reads `colors.csv` relative to the working dir, so
        # run from the repo root (where it lives) to compile past that read.
        cd(pkgdir(ShroffCelegansModels)) do
            try
                WGLMakie.activate!()
                for nerve_ring in (false, true)
                    app = App() do session::Session
                        with_theme(theme_black()) do
                            black_body(meshscatter_average(dict; session, nerve_ring))
                        end
                    end
                    mktempdir() do dir
                        export_static(joinpath(dir, "meshscatter_average.html"), app)
                    end
                end
            catch err
                @debug "MeshscatterAverage precompile workload skipped" exception = (err, catch_backtrace())
            end
        end
    end
end
end  # precompile_workload preference guard

end # module MeshscatterAverage
