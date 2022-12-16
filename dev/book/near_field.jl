
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


## this example maps the near-field E^2 around a cluster of spheres
## to compare with TERMS results

using HDF5
fid = h5open("./dev/book/map_FIELD.h5", "r")

g = fid["Near-Field"]
map_E = read(g, "map_E")
HDF5.h5_close()
x = unique(map_E[:, 2])
y = unique(map_E[:, 3])

d = (; x=map_E[:, 2], y=map_E[:, 3], z=log10.(map_E[:, 8]))
draw(data(d) * mapping(:x, :y, :z) * visual(Heatmap), axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))


using DelimitedFiles, DataFrames

positions = readdlm("./dev/book/data_field", ' ', header=false)
pos_field = DataFrame(positions, [:material, :x, :y, :z, :r])

## cluster from external data
positions = SVector.(Float64.(pos_field.x), pos_field.y, pos_field.z)
sizes = [SVector(Float64(r), r, r) for r in pos_field.r]
rotations = repeat([QuatRotation(1.0, 0.0, 0.0, 0.0)], length(positions))
materials = [String(m) for m in pos_field.material]
cl = Cluster(positions, rotations, sizes, materials, "particle")

## materials
wavelengths = [514.0]
media = Dict([("Au", epsilon_Au), ("Ag", epsilon_Ag), ("medium", x -> 1.0)])
mat = Material(wavelengths, media)

probes = SVector.(Iterators.product(x, y, zero(eltype(x))))[:]
Incidence = [RotX(0)]

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="linear"; evaluate_inside=false)

d2 = positions[.!positions.inside, :]
d2.z .= log10.(E²[.!positions.inside, 1])

# d2 = positions
# d2.z .= log10.(E²[:, 1])
# d.z .= 𝒞[.!positions.inside, 2]

draw(data(d2) * mapping(:x, :y, :z) * visual(Heatmap), axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))

terms_inside = map_E[:, 5] .!= 0
d = (; x=map_E[.!terms_inside, 2], y=map_E[.!terms_inside, 3], z=log10.(map_E[.!terms_inside, 8]))
draw(data(d) * mapping(:x, :y, :z) * visual(Heatmap), axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
