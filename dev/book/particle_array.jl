# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using DataFramesMeta
using AlgebraOfGraphics
using Makie
using Rotations
using ColorSchemes
using LaTeXStrings
home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])


## this example looks at a square array of Au nanorods in glass
## normal incidence, diffractive pitch

## materials
wavelengths = collect(400:2:1000.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## array geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl0 = cluster_single(80, 40, 40)
cl1 = cluster_array(100, 600, 80, 40, 40)
# cl1 = cluster_array(200, 600, 80, 40, 40)

## incidence: along z
Incidence = [QuatRotation(RotZYZ(0.0, 0.0, 0.0))]


disp1 = spectrum_dispersion(cl0, mat, Incidence)
disp2 = spectrum_dispersion(cl1, mat, Incidence)


d1 = dispersion_df(disp1, mat.wavelengths)
d2 = dispersion_df(disp2, mat.wavelengths)

dat1 = data(@rsubset(d1, :crosstype == "extinction"))
dat2 = data(@rsubset(d2, :crosstype == "extinction"))
m1 = dat1 * mapping(:wavelength, :value, col=:polarisation)
m2 = dat2 * mapping(:wavelength, :value, col=:polarisation)
layer1 = m1 * visual(Lines, linestyle=:dash)
layer2 = m2 * visual(Lines)
fg = draw(layer1 + layer2, facet=(; linkyaxes=:none),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))

fg

# save("figure.pdf", fg, px_per_unit=3)