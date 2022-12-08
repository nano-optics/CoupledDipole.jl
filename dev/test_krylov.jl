


push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using IterativeSolvers



using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
# using Gadfly

## variables

N_dip = 1000
N_inc = 2
N_lam = 2
N_sca = 2

wavelength = [400]
medium = Dict([("Au", CoupledDipole.epsilon_Au), ("medium", x -> 1.33)])

mat = Material(wavelength, medium)

cl = cluster_helix(N_dip, 20, 20, 30, 50, 300)

## polarizability

ii = 1;

λ = mat.wavelengths[ii]
n_medium = mat.media["medium"](λ)
kn = n_medium * 2π / λ

ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
ε = mat.media[ε_name](λ)
Alpha = alpha_particles(λ, ε, n_medium^2, cl.sizes)

ParticleRotations = map(euler_passive, cl.angles)
AlphaBlocks = map((R, A) -> R' * diagm(A) * R, ParticleRotations, Alpha)

## propagator

F = Matrix{Complex{Float64}}(I, 3 * N_dip, 3 * N_dip)
interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)

## incidence

Ejones = [SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0)]

angles = cubature_sphere(30)
N_inc = length(angles.nodes)
IncidenceRotations = map(euler_active, angles.nodes)


Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

Ein = Array{Complex{Float64}}(undef, (3 * N_dip, 2N_inc))
incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)

E = similar(Ein)
@time E = F \ Ein

Eit = copy(Ein[:, 1])
@time idrs!(Eit, F, Ein[:, 1])


test_direct = function (n)
    E = similar(Ein)
    t = @time x = F \ Ein
    t
end

test_iterative = function (n)
    for i ∈ 1:N_inc
        Eit = copy(Ein[:, i])
        idrs!(Eit, F, Ein[:, i], abstol=1e-3)
    end

end

test_direct(1)
@time test_iterative(1)
