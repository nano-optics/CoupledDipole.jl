
# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using Rotations

struct test_clust{T}

    angles::Vector{Rotation{3,T}}

end


r = rand(RotMatrix{3})
q = UnitQuaternion(r)
e = RotXYZ(0,1,2)
aa = AngleAxis(r)

test_clust([r,q,e,aa])

RotMatrix(RotZYZ(0,0,0))

RotMatrix(UnitQuaternion(1,0,0,0))

e1 = RotZYZ(0.1,0.2,0.3)
# 3×3 RotZYZ{Float64} with indices SOneTo(3)×SOneTo(3)(0.1, 0.2, 0.3):
#   0.902113  -0.383557   0.197677
#   0.387517   0.921649   0.0198338
#  -0.189796   0.0587108  0.980067


e2 = euler_active(0.1,0.2,0.3)
# 3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
#   0.902113  -0.383557   0.197677
#   0.387517   0.921649   0.0198338
#  -0.189796   0.0587108  0.980067


e1 ≈ e2

transpose(e1) ≈ euler_passive(0.1,0.2,0.3)
