
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

# rotation acting on a vector, such as an electric field

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
RotZYZ(π / 2, π / 4, π / 2) * Ex ≈ Ey
RotZYZ(π / 2, π / 4, π / 2) * Ey ≈ Ey
RotZYZ(π / 2, π / 4, π / 2) * Ez ≈ Ey


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
