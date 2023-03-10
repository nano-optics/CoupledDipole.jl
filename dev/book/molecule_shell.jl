# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFramesMeta
using DataFrames
using Rotations
using LaTeXStrings
using AlgebraOfGraphics, Makie, CairoMakie, ColorSchemes

home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
# set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])

gill(weight) = joinpath(font_folder, "GillSansNova-$(weight).otf")
set_aog_theme!(fonts=[gill("Book"), gill("Light")])

## this example looks at a spherical shell of uniaxial molecules in water
## contrasting radial and tangential configurations
## extinction spectra for varying concentrations
## orientation averaged extinction

function dye_coverage(ρ, R)
    area = 4π * R^2
    return Int(ceil(area * ρ))
end
# dye_coverage(1.2,5)

## materials
wavelength = collect(450:1:650.0)
media = Dict([("Rhodamine", alpha_lorentzmolecule), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

function model(; ρ=1, R0=2, d=0.5, medium=1.33, orientation="radial")

    R = R0 + d
    N = dye_coverage(ρ, R)
    @info "$N dipoles"
    if N > 1e4
        @warn "$N dipoles is rather a lot, are you sure? "
    end

    cl = cluster_shell(N, 0, 0, 1, R; orientation=orientation, material="Rhodamine", type="point")

    # testing
    # cl = cluster_shell(5, 1, 1, 0, R; orientation="radial", material="Rhodamine", type="point")
    res = spectrum_oa(cl, mat)
    d = oa_df(res, mat.wavelengths)
end


# model(ρ=100)
params = expand_grid(ρ=range(0.2, 1.2, step=0.2), orientation=("radial", "flat"))

all = pmap_df(params, p -> model(; p...))

## reference molecule
cl0 = cluster_single(0, 0, 1, 0, 0, 0, "Rhodamine", "point")
s = spectrum_oa(cl0, mat)
single = oa_df(s, mat.wavelengths)

d1 = data(@rsubset(all, :crosstype == "extinction", :type == "average"))
d2 = data(@rsubset(single, :crosstype == "extinction", :type == "average"))
m1 = d1 * mapping(:wavelength, :value, color=:ρ => nonnumeric, col=:orientation, row=:crosstype)
m2 = d2 * mapping(:wavelength, :value, row=:crosstype)
layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
fg = draw(layer1 + layer2, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))

fg
# save("figure.pdf", fg, px_per_unit=3)