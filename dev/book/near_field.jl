
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


## materials
wavelengths = [633.0]
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## dimer geometry
cl = cluster_chain(5, 80, 30, 30, 30, 0, 0, 0, "Au")

Incidence = [RotZ(0.0)]
## low level stuff
proto_r = cl.positions[1][1] # position type
proto_a = cl.rotations[1][1] # angle type
proto_α = 0.1 + 0.1im # dummy complex polarisability
proto_k = 2π / mat.wavelengths[1]
T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α
N_dip = length(cl.positions)
N_lam = length(mat.wavelengths)
N_inc = length(Incidence)
F = Matrix{T2}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
Ein = Array{T2}(undef, (3N_dip, 2N_inc))
E = similar(Ein)
P = similar(Ein)

Ejones = [SVector(1.0 + 0im, 0.0), SVector(0.0, 1.0 + 0im)]
# Ejones = 1.0 / √2.0 .* [
#     SVector(1im, 1.0), # Jones vector, first polar ↺
#     SVector(1.0, 1im), # Jones vector, second polar ↻
# ]

ParticleRotations = map(RotMatrix, cl.rotations) # now (active) Rotation objects
IncidenceRotations = map(RotMatrix, Incidence) # now given as quaternions
λ = mat.wavelengths[1]
n_medium = mat.media["medium"](λ)
kn = n_medium * 2π / λ
Epsilon = map(m -> mat.media[m](λ), cl.materials)
Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ)
AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)
incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)
E = F \ Ein
polarisation!(P, E, AlphaBlocks)

x = -300.0:2.0:300
# probes = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]
probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
N_pro = length(probes)


# when it was a matrix
# Esca = scattered_field(kn, cl.positions, probes, P)
# Isca = sum(reshape(abs2.(Esca), (3, N_pro * 2N_inc)), dims=1)
# plot(collect(x), log10.(Isca[1:length(x)]))

# incident field at probe locations
EinProbes = Array{T2}(undef, (3N_pro, 2N_inc))
incident_field!(EinProbes, Ejones, kn, probes, IncidenceRotations)

Esca, Bsca, Etot, Btot = local_field(kn, cl.positions, probes, P, EinProbes)

Isca = [sum(abs2.(E)) for E in Esca]
Hsca = [sum(abs2.(B)) for B in Bsca]

Itot = [sum(abs2.(E)) for E in Etot]
Htot = [sum(abs2.(B)) for B in Btot]


# plot(collect(x), log10.(Isca[1:2:length(Esca)]))
# plot(collect(x), log10.(Hsca[1:2:length(Bsca)]))
# plot(collect(x), log10.([1:2:length(Etot)]))


using HDF5
fid = h5open("./dev/book/map_Lmax1.h5", "r")
g = fid["Near-Field"]
map_E = read(g, "map_E")
map_B = read(g, "map_B")
map_C = read(g, "normalised_ldoc")
# "lambda" "x"      "y"      "z"      "scatID" "volID"  "E2avg"  "E2X"    "E2Y" 

slice = map_E[:, 3] .≈ 0.0

df = (; x=map_E[slice, 2], y=log10.(map_E[slice, 8]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=collect(x), y=log10.(Itot[1:2:length(Esca)]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)

df = (; x=map_E[slice, 2], y=log10.(map_E[slice, 9]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=collect(x), y=log10.(Itot[2:2:length(Esca)]))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)

Z₀ = 376.730313668 # free-space impedance
Y₀ = 1 / 376.730313668 # H = Y₀ E
c₀ = 299792458 # m/s

df = (; x=map_B[slice, 2], y=log10.(map_B[slice, 8]))
xy = data(df) * mapping(:x, :y)
df2 = (; x=collect(x), y=log10.(Htot[1:2:length(Esca)] / c₀ / Z₀))
xy2 = data(df2) * mapping(:x, :y)
layer = visual(Lines)
layer2 = visual(Lines, linestyle=:dash, color=:red)
draw(layer * xy + layer2 * xy2)


Esca2, Bsca2, inside = scattered_field(probes[1], kn, cl.positions, cl.sizes, ParticleRotations, P)
Esca[1:2]
Esca2

# high level


# now the near-field part 
# incident field at probe locations
Z₀ = 376.730313668 # free-space impedance
Y₀ = 1 / Z₀
T1 = eltype(P)
T2 = eltype(cl.positions[1])
Esca = [@SMatrix(zeros(T1, 3, 2N_inc)) for _ ∈ 1:N_pro]
Einc = similar(Esca)
Bsca = similar(Esca)
Etot = similar(Esca)
Btot = similar(Esca)
E² = Matrix{T2}(undef, N_pro, 2N_inc)
B² = similar(E²)
𝒞 = similar(E²)
inside = similar(E²)
for i in eachindex(probes)

    Einc[i] = incident_field(Ejones, k, probes[i], IncidenceRotations)
    Esca[i], Bsca[i], inside[i] = scattered_field(probes[i], kn, cl.positions, cl.sizes, ParticleRotations, P)
    Etot[i] = Einc[i] + Esca[i]
    Btot[i] = Y₀ * Einc[i] + Bsca[i]
    # scalar summaries
    E²[i, :] = sum(abs2.(Etot[i]), dims=1)
    B²[i, :] = sum(abs2.(Btot[i]), dims=1)
    𝒞[i, :] = imag.(sum(conj.(Etot[i]) .* Btot[i], dims=1))

end


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
