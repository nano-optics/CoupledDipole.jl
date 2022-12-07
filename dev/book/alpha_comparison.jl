# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using Rotations
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using DataFramesMeta
using VegaLite
using AlgebraOfGraphics, CairoMakie
using ColorSchemes
using LaTeXStrings
home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])


## this example compares various prescriptions of polarisability
# Mie vs Kuwata vs Majic
a
## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

a = 20.0
c = 30.0
# cl1 = cluster_dimer(200.0, a, a, c, π / 4)
cl1 = cluster_dimer(200.0, a, a, c, π / 4)

oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl1, mat; prescription="majic")
oa3 = spectrum_oa(cl1, mat; prescription="mie")

d1 = oa_df(oa1, mat.wavelengths)
d2 = oa_df(oa2, mat.wavelengths)
d3 = oa_df(oa3, mat.wavelengths)

map = mapping(:wavelength, :value, row=:crosstype, col=:type)
l1 = data(d1) * map * visual(Lines)
l2 = data(d2) * map * visual(Lines, linestyle=:dash)
l3 = data(d3) * map * visual(Lines, linestyle=:dot)
draw(l1 + l2 + l3, facet=(; linkyaxes=:none))

