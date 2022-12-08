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
home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])


## this example looks at 1 Au nanorod in water
## checking the orientation properties
# reference rod is along x
# other case is along z but rotated by pi/4, pi/2,
# with incidence along z rotated by 45deg

## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## cluster geometry

α = π / 4
β = π / 2
γ = 0.0
cl0 = cluster_single(35.0, 20.0, 20.0)
cl1 = cluster_single(20.0, 20.0, 35.0, α, β, γ)

Incidence = [RotZ(0.0)] ## incidence: along z (no rotation)
res0 = spectrum_dispersion(cl0, mat, Incidence)
d0 = dispersion_df(res0, mat.wavelengths)

Incidence = [RotZ(π / 4)] ## incidence: along z, E at 45deg
res1 = spectrum_dispersion(cl1, mat, Incidence)
d1 = dispersion_df(res1, mat.wavelengths)

map1 = mapping(:wavelength, :value, row=:polarisation, col=:crosstype)
m1 = map1 * (data(d1) * visual(Lines) +
             data(d0) * visual(Lines, linestyle=:dash, color=:red))
fg = draw(m1, facet=(; linkyaxes=:none))

fg

# save("figure.pdf", fg, px_per_unit=3)