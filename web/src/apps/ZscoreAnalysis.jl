"""
    ZscoreAnalysis

Submodule for the z-score outlier table (deployment container `zscore-analysis`,
port 9300). Like the launch-script apps, its DATA path eagerly reads PVC content
(via `launch_show_average_annotations.jl` + `zscore_analysis.jl`), so those
includes happen at **runtime** inside [`main`](@ref), never at build time.

The precompile workload exercises only the build-safe RENDER path: a synthetic
`DataFrame` with the same columns → the `:link` → `DOM.a` transform →
`Bonito.Table` → `Bonito.export_static`, caching table serialization.
"""
module ZscoreAnalysis

using WGLMakie
using Bonito
using DataFrames: DataFrame, ByRow, sort!, select!
using ShroffCelegansModels
using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference

const PORT = 9300

# Build the z-score table App from whatever `raw_zscore_analysis`/`datasets`
# resolve to at call time (defined by the runtime includes in `main`).
build_table_app() = App(; title="Shroff C. elegans z-score analysis") do
    df = raw_zscore_analysis(datasets; threshold=-Inf)
    sort!(df, :zscore, rev=true)
    select!(df, :, :link => ByRow(x->DOM.a(x; href=x)) => :link)
    return Bonito.Table(df)
end

function main(; apply_changes::Bool = true)
    # Runtime includes — these read config JSON, cell_key files, and avg_models
    # HDF5 from the PVC, so they must NOT run at precompile time.
    scripts = joinpath(pkgdir(ShroffCelegansModels), "scripts")
    Base.include(@__MODULE__, joinpath(scripts, "launch_show_average_annotations.jl"))
    Base.include(@__MODULE__, joinpath(pkgdir(ShroffCelegansModels), "src", "demo_averaging", "zscore_analysis.jl"))

    # The runtime includes define `alias_cache_unix`, `datasets`,
    # `raw_zscore_analysis`, … at a newer world age than this precompiled `main`
    # can see, so run the rest via `invokelatest`.
    Base.invokelatest() do
        @info "Activating WGLMakie"
        WGLMakie.activate!(; resize_to = :body)
        @info "Priming annotation caches"
        ShroffCelegansModels.prime_annotation_caches()
        alias_cache_unix("/nearline/shroff/")
        for group in values(datasets)
            for dataset in values(group)
                try
                    ShroffCelegansModels.load_straightened_annotations_over_time(dataset; use_myuntwist = true)
                catch err
                    @warn "Cache prime failed for dataset" path = dataset.path err
                end
            end
        end
        if apply_changes
            @info "Loading annotation changes"
            annotation_changes = ShroffCelegansModels.load_annotation_changes_cache()
            ShroffCelegansModels.update_annotations_cache(ShroffCelegansModels.annotations_cache, annotation_changes)
            @info "Loaded annotation changes"
        else
            @info "Skipping annotation changes (apply_changes = false); scoring primed cache only"
        end

        server = Server("0.0.0.0", PORT;
            proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/zscore_analysis/")
        route!(server, "/" => build_table_app())
        @info "Listening" port=PORT
        if isinteractive()
            println("Press enter to quit")
            readline()
        else
            wait(server)
        end
    end
end

# --- precompile workload (render path only) ------------------------------------

if @load_preference("precompile_workload", true)
@setup_workload begin
    df = DataFrame(
        group = ["RW10000", "RW10000"],
        embryo = [1, 2],
        annotation = ["ABprp", "ABplp"],
        timepoint = [12, 34],
        zscore = [3.2, 1.1],
        link = ["https://example.org/fix?a=1", "https://example.org/fix?a=2"],
    )
    @compile_workload begin
        try
            sort!(df, :zscore, rev=true)
            select!(df, :, :link => ByRow(x->DOM.a(x; href=x)) => :link)
            app = App(() -> Bonito.Table(df))
            mktempdir() do dir
                export_static(joinpath(dir, "zscore.html"), app)
            end
        catch err
            @debug "ZscoreAnalysis precompile workload skipped" exception = (err, catch_backtrace())
        end
    end
end
end  # precompile_workload preference guard

end # module ZscoreAnalysis
