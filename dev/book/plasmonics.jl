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
set_aog_theme!()

## this example illustrates plasmon resonances in metal nanoparticles
a
## materials
wavelength = collect(450:2:850.0)
media1 = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat1 = Material(wavelength, media1)
media2 = Dict([("Au", epsilon_Ag), ("medium", x -> 1.33)])
mat2 = Material(wavelength, media2)
media3 = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat3 = Material(wavelength, media3)

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

