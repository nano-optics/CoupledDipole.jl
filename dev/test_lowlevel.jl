# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite

## variables

N_dip = 2
N_inc = 2
N_lam = 2
N_sca = 2

wavelength = [400, 500]
medium = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])

mat = Material(wavelength, medium)


positions = [SVector(x, 0.0, 0.0) for x in (-100.0, 100.0)]
sizes = [SVector(20.0, 50.0, 20.0), SVector(20.0, 50.0, 20.0)]
rotations = repeat([QuatRotation(1.0, 0.0, 0.0, 0.0)], 2)
materials = repeat(["Au"], 2)
cl = Cluster(positions, rotations, sizes, materials, "particle")

# cl = cluster_dimer(100, 20, 20, 40, π / 4)

## polarizability

ii = 1;

λ = mat.wavelengths[ii]
n_medium = medium["medium"](λ)
kn = n_medium * 2π / λ


ε = medium["Au"](λ)

alpha_kuwata(λ, ε, n_medium^2, SVector(20, 20, 50))

map(s -> alpha_kuwata(λ, ε, n_medium^2, s), cl.sizes)

Alpha = alpha_particles(λ, ε, n_medium^2, cl.sizes)


ParticleRotations = cl.rotations
AlphaBlocks = map((R, A) -> R' * diagm(A) * R, ParticleRotations, Alpha)

## propagator

F = Matrix{Complex{Float64}}(I, 3 * N_dip, 3 * N_dip)
interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)

## incidence

Ejones = [SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0)]
Incidence = [SVector(0, 0, 0), SVector(0, π / 2, 0)]
IncidenceRotations = map(euler_active, Incidence)

euler_active(0, π / 2, 0)

Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

IncidenceRotations[2] * Evec1

Ein = Array{Complex{Float64}}(undef, (3 * N_dip, 2N_inc))
incident_field_pw!(Ein, Ejones, kn, cl.positions, IncidenceRotations)

E = similar(Ein)
E = F \ Ein

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
