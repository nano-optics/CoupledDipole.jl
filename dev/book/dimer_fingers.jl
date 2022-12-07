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


## this example looks at 2 Au nanorods in water
## contrasting positive and negative dihedral angle
## ext,scat,abs spectra as well as corresponding dichroisms
## orientation averaged

## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## dimer geometry
# cluster_dimer(d, a, b, c, ϕ=0.0, α_1=0.0, α_2=0.0, material="Au", type="particle")

cl0 = cluster_single(20.0, 20.0, 35.0)
cl1 = cluster_dimer(100.0, 20.0, 20.0, 35.0, π / 4)
cl2 = cluster_dimer(100.0, 20.0, 20.0, 35.0, -π / 4)

oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl2, mat)

d0 = oa_df(oa0, mat.wavelengths)
d1 = oa_df(oa1, mat.wavelengths)
d2 = oa_df(oa2, mat.wavelengths)

map1 = mapping(:wavelength, :value, row=:type, col=:crosstype)

m1 = map1 * (data(d1) * visual(Lines) +
             data(d2) * visual(Lines, linestyle=:dash) +
             data(d0) * visual(Lines, linestyle=:dot))
fg = draw(m1, facet=(; linkyaxes=:none))

fg

# save("figure.pdf", fg, px_per_unit=3)