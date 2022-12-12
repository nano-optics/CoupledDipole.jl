
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

x = -200.0:2.0:200
probes = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]
probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
N_pro = length(probes)

# when it was a matrix
# Esca = scattered_field(kn, cl.positions, probes, P)
# Isca = sum(reshape(abs2.(Esca), (3, N_pro * 2N_inc)), dims=1)
# plot(collect(x), log10.(Isca[1:length(x)]))


Esca = scattered_field(kn, cl.positions, probes, P)

Isca = [sum(abs2.(E)) for E in Esca]

plot(collect(x), log10.(Isca[1:2:length(Esca)]))