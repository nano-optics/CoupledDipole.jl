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

## this example compares various prescriptions of polarisability
# Mie vs Kuwata vs Majic

## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

cl0 = cluster_single(20.0, 20.0, 35.0)
cl1 = cluster_dimer(100.0, 20.0, 20.0, 35.0, π / 4)

oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl1, mat; prescription="majic")
oa3 = spectrum_oa(cl1, mat; prescription="mie")

# λ = 633.0
# e = -10.0 + 1im
# ε_m = 1.5^2
# s = SVector(1.0, 1.0, 1.0)
# CoupledDipole.alpha_majic(λ, e, ε_m, s)
# CoupledDipole.alpha_mie(λ, e, ε_m, s)
# Epsilon = [e, e]
# Sizes = [s, s]
# CoupledDipole.alpha_particles(Epsilon, Sizes, ε_m, λ; prescription="kuwata")
# CoupledDipole.alpha_particles(Epsilon, Sizes, ε_m, λ; prescription="majic")
# CoupledDipole.alpha_particles(Epsilon, Sizes, ε_m, λ; prescription="mie")



# d1 = data(@rsubset(all, :polarisation == "pol2"))
# d2 = data(@rsubset(single, :polarisation == "pol2"))

# m1 = d1 * mapping(:wavelength, :value, color=:d => nonnumeric, col=:orientation, row=:crosstype)
# m2 = d2 * mapping(:wavelength, :value, row=:crosstype)
# layer1 = m1 * visual(Lines)
# layer2 = m2 * visual(Lines, linestyle=:dash)
# draw(layer1 + layer2, facet=(; linkyaxes=:none),
#     palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))

