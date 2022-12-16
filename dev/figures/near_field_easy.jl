
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
fid = h5open("/Users/baptiste/Documents/github/book-cda/code/terms/map_easy.h5", "r")

g = fid["Near-Field"]
map_E = read(g, "map_E")
map_B = read(g, "map_B")
map_C = read(g, "normalised_ldoc")
slice = map_E[:, 3] .≈ 0.0
# "lambda" "x"      "y"      "z"      "scatID" "volID"  "E2avg"  "E2X"    "E2Y" 
HDF5.h5_close()

df = (; x=map_E[slice, 2], y=log10.(map_B[slice, 9]))
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
cl = cluster_chain(5, 80, 10, 10, 10, 0, 0, 0, "Au")
# cl = cluster_chain(5, 80, 15, 30, 30, π / 10, 0, 0, "Au")

E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="linear"; evaluate_inside=false)


slice2 = .!positions.inside
df = (; x=map_E[slice2, 2], y=log10.(map_E[slice2, 8]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=log10.(E²[:, 1]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)

df = (; x=map_E[slice2, 2], y=log10.(map_E[slice2, 9]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=log10.(E²[:, 2]))
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
df2 = (; x=positions.x, y=log10.(c₀^2 * B²[:, 1]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)


df = (; x=map_B[slice2, 2], y=(c₀^2 * map_B[slice2, 9]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=positions.x, y=(c₀^2 * B²[:, 2]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
fg = draw(layer * xy + layer2 * xy2, axis=(; xlabel="x /nm"))
# draw(layer2 * xy2)

# E², B², 𝒞, positions = map_nf(probes, cl, mat, Incidence, polarisation="circular"; evaluate_inside=false)

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


Einc, Binc, Esca, Bsca, Etot, Btot, positions = map_nf(probes[[1]], cl, mat, Incidence, polarisation="linear"; evaluate_inside=false, return_fields=true)


Etot[1]
# 1.00084+0.00277142im       0.0+0.0im
# 0.0+0.0im         0.997852-0.000721829im
# 0.0+0.0im              0.0+0.0im
# Etot
# 1.00083743370780254E+00  2.76627111331738439E-03 
# -4.11658012365589965E-19 -4.60785923306339384E-19 
# -5.00768437830872637E-20 -1.66905572781121390E-19  

# -5.55653613398821022E-19 -4.18434275943624367E-19  
# 9.97855694704689600E-01 -7.20698047446599577E-04 
# -1.53903586938544152E-21 -1.15060859056521888E-21 

# 1.00083743370780254E+00  2.76627111331738439E-03 
# -4.11658012365589965E-19 -4.60785923306339384E-19 
# -5.00768437830872637E-20 -1.66905572781121390E-19  

Etot = SVector(
    3.27068565599698915E-07 + 1im * 9.04003281775333817E-10,
    9.97855694704636309E-01 + 1im * -7.20698047446561088E-04,
    -1.53905223424234491E-21 + 1im * -1.15066313445454197E-21)

Esca[1]
# 0.000839231+0.00277142im          0.0+0.0im
# 0.0+0.0im         -0.00214833-0.000721829im
# 0.0+0.0im                 0.0+0.0im
# Esca
# 8.37433707802585914E-04  2.76627111331738439E-03 
# -4.11658012365589965E-19 -4.60785923306339384E-19 
# -5.00768437830872637E-20 -1.66905572781121390E-19

# -5.55653613398821022E-19 -4.18434275943624367E-19 
# -2.14430529531045751E-03 -7.20698047446599577E-04 
# -1.53903586938544152E-21 -1.15060859056521888E-21

Btot[1]
# 0.0+0.0im   -4.4364e-9+0.0im
# 4.4364e-9+0.0im          0.0+0.0im
#       0.0+0.0im  1.11867e-11+8.6553e-13im

# 1.16779024513335779E-43  3.16019347520232125E-44  
# 4.43640246613542267E-09  1.04639635433630616E-29  
# 1.50328100036902131E-27 -8.53934769440390897E-28  

# -4.43640246613542267E-09  5.29654502546095839E-29 
# -7.93270520758074862E-30 -4.32507174903303460E-30  
# 1.11659763155243181E-11  8.64991445557789858E-13  

# 1.16779024513335779E-43  3.16019347520232125E-44  
# 4.43640246613542267E-09  1.04639635433630616E-29  
# 1.50328100036902131E-27 -8.53934769440390897E-28  



Btot = SVector(
    -4.43640246613518609E-09 + 1im * 5.29654502546066580E-29,
    1.44979368492225930E-15 + 1im * -4.32506832946291867E-30,
    1.11659763155237219E-11 + 1im * 8.64991445557744218E-13
)



Bsca[1]
# 0.0+0.0im          0.0+0.0im
# 0.0+0.0im          0.0+0.0im
# 0.0+0.0im  1.11867e-11+8.6553e-13im

# 1.16779024513335779E-43  3.16019347520232125E-44  
# 7.40968438194652487E-30  1.04639635433630616E-29  
# 1.50328100036902131E-27 -8.53934769440390897E-28

# 6.83718857708100834E-28  5.29654502546095839E-29 
# -7.93270520758074862E-30 -4.32507174903303460E-30  
# 1.11659763155243181E-11  8.64991445557789858E-13


# Etot1
# 1.00083743370780254E+00  2.76627111331738439E-03 
# -4.11658012365589965E-19 -4.60785923306339384E-19 
# -5.00768437830872637E-20 -1.66905572781121390E-19  

# 8.37433707802585914E-04  2.76627111331738439E-03 
# -4.11658012365589965E-19 -4.60785923306339384E-19 
# -5.00768437830872637E-20 -1.66905572781121390E-19  
# Etot2
# 3.27068565599698915E-07  9.04003281775333817E-10 
# 9.97855694704636309E-01 -7.20698047446561088E-04 
# -1.53905223424234491E-21 -1.15066313445454197E-21  

# 2.73669061560564628E-10  9.04003281775333817E-10 
# -2.14430529531034302E-03 -7.20698047446561088E-04 
# -1.53905223424234491E-21 -1.15066313445454197E-21


# C
# -1.13952372562485136E-09