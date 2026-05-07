using Dates

function (@main)(args)
    date = Dates.format(Dates.today(), "yyyy_mm_dd")

    pod = strip(read(`oc get pod -n shroff-data -l app=shroff-data -o $("jsonpath={.items[0].metadata.name}")`, String))

    dest = joinpath(@__DIR__, "..", "annotation_changes_$(date).h5")

    run(`oc cp -n shroff-data -c fix-ap-axis $(pod):/data/annotations/annotation_changes.h5 $(dest)`)

    println("Downloaded to $dest")
end
