
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


## this example maps the near-field E^2 around a nanolens of nanospheres
## to compare with TERMS results

using HDF5
# fid = h5open("./dev/book/map_Lmax1.h5", "r")
fid = h5open("$home/Documents/github/book-cda/code/terms/map_nanolens_dipole.h5", "r")

g = fid["Near-Field"]
map_E = read(g, "map_E")
HDF5.h5_close()
# "lambda" "x"      "y"      "z"      "scatID" "volID"  "E2avg"  "E2X"    "E2Y" 

df0 = (; x=map_E[:, 2], y=map_E[:, 3], z=log10.(map_E[:, 7]))

using HDF5
# fid = h5open("./dev/book/map_Lmax1.h5", "r")
fid = h5open("$home/Documents/github/book-cda/code/terms/map_nanolens.h5", "r")

g = fid["Near-Field"]
map_E = read(g, "map_E")
HDF5.h5_close()
# "lambda" "x"      "y"      "z"      "scatID" "volID"  "E2avg"  "E2X"    "E2Y" 

df1 = (; x=map_E[:, 2], y=map_E[:, 3], z=log10.(map_E[:, 7]))

layer = visual(Heatmap)
map = mapping(:x, :y, :z)

f1 = data(df1) * layer * map
f0 = data(df0) * layer * map
draw(f1, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
# draw()

## materials
wavelengths = [369.0]
media = Dict([("Au", epsilon_Ag), ("medium", x -> 1.0)])
mat = Material(wavelengths, media)

Incidence = [RotZ(0.0)]
# x = range(-50.0, 50, 500)
# y = range(0.0, 100, 500)
# probes = SVector.(Iterators.product(x, y, zero(eltype(x))))[:]
probes = SVector.(df0.x, df0.y, 0.0)

positions = [SVector(0.0, 0, 0), SVector(0.0, 69, 0), SVector(0.0, 92, 0)]
r1 = 45
r2 = 15
r3 = 5
sizes = [SVector(r1, r1, r1), SVector(r2, r2, r2), SVector(r3, r3, r3)]
rotations = repeat([QuatRotation(1.0, 0.0, 0.0, 0.0)], 3)
materials = repeat(["Au"], 3)
cl = Cluster(positions, rotations, sizes, materials, "particle")

# cl = cluster_chain(5, 80, 15, 30, 30, π / 10, 0, 0, "Au")

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="linear"; evaluate_inside=false)


df2 = (; x=positions.x, y=positions.y, z=log10.(E²[:, 2]))
f2 = data(df2) * layer * map
draw(f2, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()); figure=(resolution=(600, 400),))


using Makie
resolution = (800, 400)
fig = Figure(; resolution)
draw!(fig[1, 1], f1, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
draw!(fig[1, 2], f0, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
draw!(fig[1, 3], f2, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
fig

## errors

# log of relative error
df4 = (; x=positions.x, y=positions.y, z=log10.(abs.(10 .^ df0.z - 10 .^ df2.z) ./ abs.(10 .^ df0.z)))
df5 = (; x=positions.x, y=positions.y, z=log10.(abs.(10 .^ df1.z - 10 .^ df2.z) ./ abs.(10 .^ df1.z)))
f4 = data(df4) * layer * map
f5 = data(df5) * layer * map


draw(f4, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
draw(f5, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))

resolution = (800, 400)
fig = Figure(; resolution)
draw!(fig[1, 1], f4, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
draw!(fig[1, 2], f5, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
fig

