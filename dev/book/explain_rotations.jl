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
using AlgebraOfGraphics, Makie, CairoMakie
home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
# set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])

gill(weight) = joinpath(font_folder, "GillSansNova-$(weight).otf")
set_aog_theme!(fonts=[gill("Book"), gill("Light")])

## this example illustrates rotations
# as implemented in Rotations.jl
# but also with standard Euler rotations


# Active Euler rotation in ZYZ convention
α = π / 3;
β = π / 4;
γ = π / 5;
R1 = RotZYZ(α, β, γ)
RotMatrix(R1)

R2 = euler_active(α, β, γ)
R1 ≈ R2

R2p = euler_passive(α, β, γ)
R2p ≈ R1 # false
R2p ≈ inv(R1)
R2p ≈ transpose(R1)
R2p ≈ R1'
# , such as an electric field

# rotation acting on a vector

E1 = SVector(1.0, 0.0, 0.0)
E2 = 1 / sqrt(2) .* SVector(0.0 + 1im, 1.0, 0.0)
R1 * E1
R1 * E2

# OK, but how do we know these are correct?
# let's try some obvious cases

Ex = SVector(1.0, 0.0, 0.0)
Ey = SVector(0.0, 1.0, 0.0)
Ez = SVector(0.0, 0.0, 1.0)
RotZYZ(π / 2, 0.0, 0.0) * Ex ≈ Ey
RotZYZ(π / 2, 0.0, 0.0) * Ey ≈ -Ex
RotZYZ(π / 2, 0.0, 0.0) * Ez ≈ Ez
RotZYZ(0.0, π / 2, 0.0) * Ex ≈ -Ez
RotZYZ(0.0, π / 2, 0.0) * Ey ≈ Ey
RotZYZ(0.0, π / 2, 0.0) * Ez ≈ Ex
# now combine 3 rotations
RotZYZ(π / 2, π / 2, π / 4) * Ex ≈ 1 / √2 .* (-Ex - Ez) # Ex -> Ey -> -Ez -> 1/√2 (-Ex-Ez)
RotZYZ(π / 2, π / 2, π / 4) * Ey ≈ 1 / √2 .* (-Ex + Ez)
RotZYZ(π / 2, π / 2, π / 4) * Ez ≈ Ey

# sounds reasonable, now how do we apply these?

# for incident and scattered directions, we want to actively rotate vectors in the common lab frame
# so we use directly active rotations 

# in spectrum_dispersion
# ParticleRotations = map(RotMatrix, cl.rotations)
# IncidenceRotations = map(RotMatrix, Incidence) 
# ScatteringVectors = map(euler_unitvector, quad_sca.nodes)

# incident_field!
# rotating the electric field 
# Rm = IncidenceRotations[jj]
# Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
# E1_r = Rm * Evec1 # active rotation
# k_hat = kn * Rm[:, 3] # Rm applied to kx=(0,0,1) gives third column, then scaled by scalar kn

# for particles, however, we want to transform between local frame and lab frame
# so a passive transformation is used (describing the same object in two different coordinate systems)

# spectrum_dispersion, spectrum_oa
# rotating the polarisability tensor from its particle frame to the lab frame
# if R is passive:
# AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
# if R was active:
# AlphaBlocks = map((R, A) -> R * (diagm(A) * R'), ParticleRotations, Alpha)

# clusters
# rotations expresses the orientation of the particle 
# in the lab frame, and is the inverse of the active rotation to rotate the particle
# because we want to express the passive rotation associated with the frame transformation

# example, in cluster_single:
# rotations = [inv(QuatRotation(Rotations.RotZYZ(α, β, γ)))]

## testing with a single particle
# cf single_rod.jl






q = QuatRotation(r)
e = RotXYZ(0, 1, 2)
aa = AngleAxis(r)



test_clust([r, q, e, aa])


RotMatrix(QuatRotation(1, 0, 0, 0))

e1 = RotZYZ(0.1, 0.2, 0.3)
# 3×3 RotZYZ{Float64} with indices SOneTo(3)×SOneTo(3)(0.1, 0.2, 0.3):
#   0.902113  -0.383557   0.197677
#   0.387517   0.921649   0.0198338
#  -0.189796   0.0587108  0.980067


e2 = euler_active(0.1, 0.2, 0.3)
# 3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
#   0.902113  -0.383557   0.197677
#   0.387517   0.921649   0.0198338
#  -0.189796   0.0587108  0.980067


e1 ≈ e2

transpose(e1) ≈ euler_passive(0.1, 0.2, 0.3)
