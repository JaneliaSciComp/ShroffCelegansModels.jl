using WGLMakie
using Bonito
using Sockets


if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

include("../../scripts/launch_show_average_annotations.jl")
include("../../src/demo_averaging/zscore_analysis.jl")

function web_zscore_analysis(datasets = datasets)
    @info "Loading annotation changes"
    annotation_changes = ShroffCelegansModels.load_annotation_changes_cache()
    ShroffCelegansModels.update_annotations_cache(ShroffCelegansModels.annotations_cache, annotation_changes);
    @info "Loaded annotation changes"
    shroff_data_ip = "0.0.0.0"
    server = Server(
        shroff_data_ip, 9300;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/zscore_analysis/"
    )
    route!(server, "/" => App(; title="Shroff C. elegans z-score analysis") do
        df = raw_zscore_analysis(datasets)
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