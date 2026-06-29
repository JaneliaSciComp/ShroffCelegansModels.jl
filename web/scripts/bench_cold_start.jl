# Cold-start benchmark for the web apps' render paths.
#
# Run in a FRESH Julia process so the timing reflects what a container sees on
# startup (the first render is where Makie → WGLMakie → Bonito compilation hits):
#
#     julia --project=web web/scripts/bench_cold_start.jl [app]
#
# where [app] is one of: meshscatter (default), modified, zscore, common.
#
# With the PrecompileTools workload baked into the image the first render is a
# small fraction of its cost without the cache. A/B by toggling the
# `precompile_workload` Preference and re-precompiling (cache cleared between).
using ShroffCelegansModelsWebInterface
using ShroffCelegansModelsWebInterface: MeshscatterAverage, ModifiedTimes, ZscoreAnalysis, CommonScenes
using Bonito: export_static, App, Table, DOM
using DataFrames: DataFrame, ByRow, sort!, select!
using WGLMakie

# Returns a thunk that builds the App for the requested app's render path,
# using the same synthetic inputs the corresponding @compile_workload uses.
function app_builder(which::AbstractString)
    if which == "meshscatter"
        dict = MeshscatterAverage._synthetic_annotation_dict()
        return () -> (first(MeshscatterAverage.build_apps(dict)))
    elseif which == "modified"
        kinds = ModifiedTimes._synthetic_kinds()
        return () -> App(() -> ModifiedTimes.render_kinds(kinds, 1.7e9))
    elseif which == "zscore"
        df = DataFrame(group=["g","g"], embryo=[1,2], annotation=["a","b"],
                       timepoint=[12,34], zscore=[3.2,1.1],
                       link=["https://e.org/1","https://e.org/2"])
        sort!(df, :zscore, rev=true)
        select!(df, :, :link => ByRow(x->DOM.a(x; href=x)) => :link)
        return () -> App(() -> Table(df))
    elseif which == "common"
        return () -> App(() -> CommonScenes._representative_figure())
    else
        error("unknown app '$which'")
    end
end

function main(which = isempty(ARGS) ? "meshscatter" : ARGS[1])
    build = app_builder(which)
    WGLMakie.activate!()
    t_render = @elapsed begin
        app = build()
        mktempdir() do dir
            export_static(joinpath(dir, "out.html"), app)
        end
    end
    @info "time-to-first-render" app = which seconds = t_render
end

main()
