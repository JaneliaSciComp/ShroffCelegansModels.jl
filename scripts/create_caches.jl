for (key,dataset) in datasets
    @info "showing" key
    try
        fig = show_average_annotations(avg_models, dataset; use_myuntwist=true);
    catch err
    end
    fig = nothing
    GC.gc()
    GC.gc()
    GC.gc()
    GC.gc()
    sleep(2)
end