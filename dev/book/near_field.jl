
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


## this example maps the near-field E^2 around a chain of nanospheres
## to compare with TERMS results

using HDF5
# fid = h5open("./dev/book/map_Lmax1.h5", "r")
fid = h5open("/Users/baptiste/Documents/github/book-cda/code/terms/map_Lmax1.h5", "r")

g = fid["Near-Field"]
map_E = read(g, "map_E")
map_B = read(g, "map_B")
map_C = read(g, "normalised_ldoc")
slice = map_E[:, 3] .≈ 0.0
# "lambda" "x"      "y"      "z"      "scatID" "volID"  "E2avg"  "E2X"    "E2Y" 

df = (; x=map_E[slice, 2], y=log10.(map_B[slice, 8]))
xy = data(df) * mapping(:x, :y)
layer = visual(Lines)
draw(layer * xy)


## materials
wavelengths = [633.0]
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## line cut
Incidence = [RotZ(0.0)]
x = -300.0:1.0:300
y = -100.0:1.0:100
probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
cl = cluster_chain(5, 80, 30, 30, 30, 0, 0, 0, "Au")
# cl = cluster_chain(5, 80, 15, 30, 30, π / 10, 0, 0, "Au")

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="circular"; evaluate_inside=false)


slice2 = .!positions.inside
df = (; x=map_E[slice, 2], y=log10.(map_E[slice, 8]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=log10.(E²[:, 2]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)

df = (; x=map_E[slice2, 2], y=log10.(map_E[slice2, 9]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=log10.(E²[:, 1]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)

Z₀ = 376.730313668 # free-space impedance
Y₀ = 1 / 376.730313668 # H = Y₀ E
c₀ = 299792458 # m/s

slice2 = .!positions.inside
df = (; x=map_B[slice2, 2], y=log10.(c₀^2 * map_B[slice2, 8]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=log10.(c₀^2 * B²[:, 2]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)


df = (; x=map_B[slice2, 2], y=(c₀^2 * map_B[slice2, 9]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=(c₀^2 * B²[:, 1]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
fg = draw(layer * xy + layer2 * xy2, axis=(; xlabel="x /nm"))
# draw(layer2 * xy2)

probes2 = probes[[1]]

Einc, Binc, Esca, Bsca, Etot, Btot, positions = map_nf(probes2, cl, mat, Incidence, polarisation="linear"; evaluate_inside=false, return_fields=true)
Etot[1]
Einc[1]
Esca[1]
Binc[1]
Bsca[1]
Btot[1]


# Etot
# 9.50215901807303531E-01  1.18776253408519855E-01 
# -7.80625564189563192E-18 -2.94902990916057206E-17  
# 1.22044526811645668E-03  1.68799142884323752E-04 
# vs 
# 0.945584+0.121082im       0.0+0.0im
# 0.0+0.0im       0.950494-0.0221498im
# 0.0+0.0im            0.0+0.0im
# Esca
# -4.97840981926966630E-02  1.18776253408519855E-01 
# -7.80625564189563192E-18 -2.94902990916057206E-17  
# 1.22044526811645668E-03  1.68799142884323752E-04
# vs
# -0.0544159+0.121082im         0.0+0.0im
# 0.0+0.0im       -0.0495062-0.0221498im
# 0.0+0.0im              0.0+0.0im

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="circular"; evaluate_inside=false)

df = (; x=map_C[slice2, 2], y=(map_C[slice2, 6]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=(𝒞[:, 1]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
fg = draw(layer * xy + layer2 * xy2, axis=(; xlabel="x /nm"))


df = (; x=map_C[slice2, 2], y=(map_C[slice2, 7]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=(𝒞[:, 2]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
fg = draw(layer * xy + layer2 * xy2, axis=(; xlabel="x /nm"))


# high level


# now the near-field part 
# incident field at probe locations

x = -300.0:1.0:300
y = -100.0:1.0:100
probes = SVector.(Iterators.product(x, y, zero(eltype(x))))[:]
Incidence = [RotX(0 * π / 2)]

cl = cluster_chain(5, 80, 30, 30, 30, 0, 0, 0, "Au")
# cl = cluster_chain(5, 80, 15, 30, 30, π / 10, 0, 0, "Au")

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="circular")

d = positions[.!positions.inside, :]
d.z .= 𝒞[.!positions.inside, 1]

draw(data(d) * mapping(:x, :y, :z) * visual(Heatmap), axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()))
