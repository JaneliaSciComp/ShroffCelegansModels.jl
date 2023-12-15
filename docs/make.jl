using ShroffCelegansModels
using Documenter

DocMeta.setdocmeta!(ShroffCelegansModels, :DocTestSetup, :(using ShroffCelegansModels); recursive=true)

makedocs(;
    modules=[ShroffCelegansModels],
    authors="Mark Kittisopikul <kittisopikulm@janelia.hhmi.org> and contributors",
    repo="https://github.com/Mark Kittisopikul/ShroffCelegansModels.jl/blob/{commit}{path}#{line}",
    sitename="ShroffCelegansModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Mark Kittisopikul.github.io/ShroffCelegansModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Mark Kittisopikul/ShroffCelegansModels.jl",
    devbranch="main",
)
