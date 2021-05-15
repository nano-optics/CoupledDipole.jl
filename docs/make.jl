using CoupledDipole
using Documenter

DocMeta.setdocmeta!(CoupledDipole, :DocTestSetup, :(using CoupledDipole); recursive=true)

makedocs(;
    modules=[CoupledDipole],
    authors="baptiste",
    repo="https://github.com/nano-optics/CoupledDipole.jl/blob/{commit}{path}#{line}",
    sitename="CoupledDipole.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nano-optics.github.io/CoupledDipole.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nano-optics/CoupledDipole.jl.git",
)
