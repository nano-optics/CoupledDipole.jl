using CoupledDipole
using Documenter
using StaticArrays

DocMeta.setdocmeta!(CoupledDipole, :DocTestSetup, :(using CoupledDipole, StaticArrays); recursive=true)

makedocs(;
    modules=[CoupledDipole],
    authors="baptiste",
    repo="https://github.com/nano-optics/CoupledDipole.jl/blob/{commit}{path}#{line}",
    sitename="CoupledDipole.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nano-optics.github.io/CoupledDipole.jl",
        assets=[asset("assets/three.min.js", class=:js, islocal=true),
        asset("assets/Detector.js", class=:js, islocal=true),
        asset("assets/TrackballControls.js", class=:js, islocal=true)],
    ),
    pages=[
        "Home" => "index.md",
        "gettingstarted.md",
        "Introduction" => [
            "theory.md"
        ],
        "Examples" => [
            "dimer.md",
            "helix.md",
            "array.md",
            "dyes.md"
        ],
        "Implementation" => [
            "conventions.md",
            "high_level.md",
            "clusters.md",
            "materials.md",
            "system.md",
            "cross_sections.md",
            "utils.md",
            "visualise.md"
        ],
    ],
)


deploydocs(;
    repo="github.com/nano-optics/CoupledDipole.jl.git",
    devbranch = "main"
)
