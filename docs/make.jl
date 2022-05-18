using FiniteVolumeFluid
using Documenter

DocMeta.setdocmeta!(FiniteVolumeFluid, :DocTestSetup, :(using FiniteVolumeFluid); recursive=true)

makedocs(;
    modules=[FiniteVolumeFluid],
    authors="iamzhtr@hotmail.com",
    repo="https://github.com/HetaoZ/FiniteVolumeFluid.jl/blob/{commit}{path}#{line}",
    sitename="FiniteVolumeFluid.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/FiniteVolumeFluid.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/FiniteVolumeFluid.jl",
    devbranch="main",
)
