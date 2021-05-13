# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
# using Gadfly

## variables

N_dip = 2
N_inc = 2
N_lam = 2
N_sca = 2

wavelength = [400, 500]
medium = Dict([("Au", CoupledDipole.epsilon_Au), ("medium", x -> 1.33)])

mat = Material(wavelength, medium)

cl = cluster_dimer(100, 20, 20, 40, π/4)

## polarizability

ii=1;

λ = mat.wavelength[ii]
n_medium = mat.medium["medium"](λ)
kn = n_medium * 2π / λ

ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
ε = mat.medium[ε_name](λ)
Alpha = alpha_spheroids(λ, ε, n_medium^2, cl.sizes)

ParticleRotations = map(euler_passive, cl.angles)
AlphaBlocks = map((R, A) -> R' * diagm(A) * R, ParticleRotations, Alpha)

## propagator

F = Matrix{Complex{Float64}}(I, 3*N_dip, 3*N_dip)
propagator_freespace_labframe!(F,kn,cl.positions,AlphaBlocks)

## incidence

Ejones = [SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0)]
Incidence = [SVector(0,0,0),SVector(0,π/2,0)]
IncidenceRotations = map(euler_active, Incidence)

euler_active(0,π/2,0)

Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

IncidenceRotations[2] * Evec1

Ein = Array{Complex{Float64}}(undef, (3*N_dip, 2N_inc))
incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)

E = similar(Ein)
E = F\Ein

P = similar(Ein)
polarisation!(P, E, AlphaBlocks)


## cross sections

Cext = Array{Float64}(undef, (2N_inc))
Cabs = similar(Cext)
Csca = similar(Cext)

extinction!(Cext, kn, P, Ein)
absorption!(Cabs, kn, P, E)

q = cubature_sphere(36)
scattering!(Csca, cl.positions, q.nodes, q.weights, kn, P)
