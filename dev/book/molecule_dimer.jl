# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using Rotations
using DataFramesMeta
using ColorSchemes
set_aog_theme!()

## this example looks at 2 uniaxial molecules in water
## contrasting head-t-t and side-b-s configurations
## extinction spectra for varying distances
## fixed orientation, same as Au NP example but different scale


## materials
wavelength = collect(450:1:650.0)
media = Dict([("Rhodamine", alpha_bare), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

function model(; d=100, orientation="head-to-tail")
    if orientation == "head-to-tail"
        # dimer axis along y
        positions = [SVector(0.0, y, 0.0) for y in (-d / 2.0, d / 2.0)]
    elseif orientation == "side-by-side"
        # dimer axis along x
        positions = [SVector(x, 0.0, 0.0) for x in (-d / 2.0, d / 2.0)]
    end

    sizes = repeat([SVector(0, 1.0, 0.0)], 2)
    rotations = repeat([QuatRotation(1.0, 0.0, 0.0, 0.0)], 2)
    materials = repeat(["Rhodamine"], 2)
    cl = Cluster(positions, rotations, sizes, materials, "point")

    Incidence = [RotZ(0.0)] ## incidence: along z (no rotation)
    res = spectrum_dispersion(cl, mat, Incidence)
    d = dispersion_df(res, mat.wavelengths)
end


td = [collect(range(0.7, 1.3, step=0.2)); 5]
params = expand_grid(d=td, orientation=("head-to-tail", "side-by-side"))

all = pmap_df(params, p -> model(; p...))

## reference molecule
cl0 = cluster_single(0, 1, 0, 0, 0, 0, "Rhodamine", "point")
s = spectrum_dispersion(cl0, mat, [QuatRotation(RotZ(0.0))])
single = dispersion_df(s, mat.wavelengths)


d1 = data(@rsubset(all, :polarisation == "pol2"))
d2 = data(@rsubset(single, :polarisation == "pol2"))
m1 = d1 * mapping(:wavelength, :value, color=:d => nonnumeric, col=:orientation, row=:crosstype)
m2 = d2 * mapping(:wavelength, :value, row=:crosstype)
layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
fg = draw(layer1 + layer2, facet=(; linkyaxes=:none),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))

fg

# save("figure.pdf", fg, px_per_unit=3)