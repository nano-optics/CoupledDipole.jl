
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



## materials
wavelengths = [633.0]
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])

media = Dict([("dummy", x -> 1.33^2), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## line cut
Incidence = [RotZ(0.0)]
x = -300.0:1.0:300
φ = range(0, 2π, 36)
θ = range(0, π, 18)
directions = SVector.(Iterators.product(φ, θ, 0.0))[:]

directions = SVector.(Iterators.product(0.0, θ, 0.0))[:]
cl = cluster_chain(5, 80, 30, 30, 30, 0, 0, 0, "dummy")
source = SVector(0, 10.0, 0.0) # not at origin...

E² = scattering_pattern(directions,
    cl::Cluster,
    mat::Material,
    source;
    prescription="kuwata")


df = (; x=θ, y1=log10.(E²[:, 1]), y2=log10.(E²[:, 2]), y3=log10.(E²[:, 3]))
df = (; x=θ, y1=(E²[:, 1]), y2=(E²[:, 2]), y3=(E²[:, 3]))
xy = data(df) * mapping(:x, :y1)
xy2 = data(df) * mapping(:x, :y2)
xy3 = data(df) * mapping(:x, :y3)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
layer3 = visual(Lines, linestyle=:dot, color=:blue)
draw(layer * xy + layer2 * xy2 + layer3 * xy3)
