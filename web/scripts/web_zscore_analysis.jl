using WGLMakie
using Bonito
using Sockets


if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

include("../../scripts/launch_show_average_annotations.jl")
include("../../src/demo_averaging/zscore_analysis.jl")

"""
    web_zscore_analysis(datasets = datasets; apply_changes = true)

Serve the z-score outlier table.

The score is computed from `annotations_cache` (via `raw_annotation_dict` →
`load_straightened_annotations_over_time`). For edits to be reflected, the cache
must be (1) primed for every dataset so `update_annotations_cache` has keys to
update, then (2) updated with the live annotation changes. Without priming, the
update applies to an empty cache (every change is skipped) and the table scores
raw, unedited positions.

`apply_changes = true` (default) primes from the recompute output (which already
has the last run's edits baked in) and then re-applies the live
`annotation_changes.h5` so edits made since the last recompute are included.
Set `apply_changes = false` to score only the recompute pipeline's end-of-run
state (the primed cache) without layering the most recent changes on top.
"""
function web_zscore_analysis(datasets = datasets; apply_changes::Bool = true)
    # Prime annotations_cache from the freshest available HDF5 (recompute PVC
    # output, else baked-in snapshot), then map any legacy Windows-rooted keys
    # to the Linux dataset.path. Mirrors web_show_average_annotations.jl.
    @info "Priming annotation caches"
    ShroffCelegansModels.prime_annotation_caches()
    alias_cache_unix("/nearline/shroff/")
    # Ensure every dataset is present as a cache key before applying changes —
    # update_annotations_cache only mutates existing entries (recompute phase 2).
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
        ShroffCelegansModels.update_annotations_cache(ShroffCelegansModels.annotations_cache, annotation_changes);
        @info "Loaded annotation changes"
    else
        @info "Skipping annotation changes (apply_changes = false); scoring primed cache only"
    end
    shroff_data_ip = "0.0.0.0"
    server = Server(
        shroff_data_ip, 9300;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/zscore_analysis/"
    )
    route!(server, "/" => App(; title="Shroff C. elegans z-score analysis") do
        df = raw_zscore_analysis(datasets; threshold=-Inf)
        sort!(df, :zscore, rev=true)
        select!(df, :, :link => ByRow(x->DOM.a(x; href=x)) => :link)
        table = Bonito.Table(df)
        return table
    end)
    return server
end

function web_main()
    WGLMakie.activate!(; resize_to = :body)
    server = web_zscore_analysis()
end

println(@__FILE__)
println(abspath(PROGRAM_FILE))

if run_web_main
    wait(web_main())
    println("Press enter to quit")
    readline()
end
