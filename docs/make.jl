using MeanFieldGraph
using Documenter

DocMeta.setdocmeta!(MeanFieldGraph, :DocTestSetup, :(using MeanFieldGraph); recursive=true)

makedocs(;
    modules=[MeanFieldGraph],
    authors="Julien Chevallier",
    sitename="MeanFieldGraph.jl",
    format=Documenter.HTML(;
        canonical="https://jucheval.github.io/MeanFieldGraph.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jucheval/MeanFieldGraph.jl",
    devbranch="main",
)
