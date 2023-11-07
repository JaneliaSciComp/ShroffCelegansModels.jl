using ShroffCelegansModel
using Documenter

DocMeta.setdocmeta!(ShroffCelegansModel, :DocTestSetup, :(using ShroffCelegansModel); recursive=true)

makedocs(;
    modules=[ShroffCelegansModel],
    authors="Mark Kittisopikul <kittisopikulm@janelia.hhmi.org> and contributors",
    repo="https://github.com/Mark Kittisopikul/ShroffCelegansModel.jl/blob/{commit}{path}#{line}",
    sitename="ShroffCelegansModel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Mark Kittisopikul.github.io/ShroffCelegansModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Mark Kittisopikul/ShroffCelegansModel.jl",
    devbranch="main",
)
