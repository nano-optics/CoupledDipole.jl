using Rotations


function dimer_rotations(ϕ = 0.0, α_1 = 0.0, α_2 = 0.0)

    q1 = QuatRotation(cos(α_1/2), sin(α_1/2), 0, 0) # rotation α_1 about x
    q2 = QuatRotation(cos(α_2/2), sin(α_2/2), 0, 0) # rotation α_2 about x
    q3 = QuatRotation(cos(ϕ/2), 0, sin(ϕ/2), 0) # rotation ϕ about y
    # rotate particle 1 by q1 only (stays in yz plane)
    # rotate particle 2 by q2, then q3 but in original frame so order swapped
    return [q1, q3*q2]
end

dimer_rotations(0.0, 0.0, 0.0)
# identities, OK

a = dimer_rotations(pi/2, 0.0, 0.0)
# 2-element Vector{QuatRotation{Float64}}:
#  [1.0 0.0 0.0;
#   0.0 1.0 0.0;
#   0.0 0.0 1.0]
#  [0.0 0.0 1.0;
#   0.0 1.0 0.0;
#   -1.0 0.0 0.0]

using StaticArrays
a[2] * SVector(1,0,0) ≈ SVector(0,0,-1) # x becomes -z
a[2] * SVector(0,1,0) ≈ SVector(0,1,0) # y stays y
a[2] * SVector(0,0,1) ≈ SVector(1,0,0) # z becomes x

b = dimer_rotations(pi/4, 0.0, 0.0)
b[2] * SVector(1,0,0) ≈ sqrt(2)/2 .* SVector(1,0,-1)
b[2] * SVector(0,1,0) ≈ SVector(0,1,0)
b[2] * SVector(0,0,1) ≈ sqrt(2)/2 .* SVector(1,0,1)

c = dimer_rotations(pi/4, 0.0, pi/2 - pi/100)
# first rotation around x by almost 90 degrees
# then along original y by pi/4
c[2] * SVector(1,0,0) # x stays x, then √2/2(1,0,-1)
c[2] * SVector(0,1,0) # y becomes almost z, then √2/2(1,0,1)
c[2] * SVector(0,0,1) # z becomes almost -y, then almost stays there
