
## cluster definitions


"""
    Cluster(positions, angles, sizes, material)

Particle cluster specification.

"""
struct Cluster{T1,T2,T3}

    "positions::Vector{SVector{3,T1}}"
    positions::Vector{SVector{3,T1}}

    "rotations::Vector{UnitQuaternion{T2}}"
    rotations::Vector{UnitQuaternion{T2}}

    "sizes::Vector{SVector{3,T3}}"
    sizes::Vector{SVector{3,T3}}

    "material::Vector{String}"
    material::Vector{String}

    "type::Vector{String}"
    type::Vector{String}
end



@doc raw"""
    cluster_single(a::T, b::T, c::T, α::T = 0.0, β::T = 0.0, γ::T = 0.0) where T <: Real

Particle cluster consisting of a single particle at the origin
- `a,b,c`: semi-axis along x,y,z
- `α,β,γ`: Euler angles
- `material`: String referencing the material of the particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_single(1.0,2.0,3.0)
```

"""
function cluster_single(a, b, c, α = 0.0, β = 0.0, γ = 0.0, material = "Au", type="particle")
    sizes = [SVector(a, b, c)]
    positions = [SVector(0.0, 0.0, 0.0)]
    # input parameters are Euler angles
    rotations = [UnitQuaternion(RotZYZ(α,β,γ))]
    Cluster(positions, rotations, sizes, [material], [type])
end


"""
    cluster_dimer(d::T, a::T, b::T, c::T, ϕ::T = 0.0, α_1::T = 0.0, α_2::T = 0.0) where T <: Real

Particle cluster consisting of 2 identical particles separated along y
- `a,b,c`: semi-axes along x,y,z
- `ϕ`: dihedral angle between both particles seen along the y-axis
- `α_1,α_2`: angle of each particle with the y axis
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_dimer(80, 10, 10, 20)
```

"""
function cluster_dimer(d, a, b, c, ϕ = 0.0, α_1 = 0.0, α_2 = 0.0, material = "Au", type="particle")
    sizes = [SVector(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector(0.0, y, 0.0) for y in (-d/2, d/2)]
    q1 = UnitQuaternion(cos(α_1/2), sin(α_1/2), 0, 0) # rotation α_1 about x
    q2 = UnitQuaternion(cos(α_2/2), sin(α_2/2), 0, 0) # rotation α_2 about x
    q3 = UnitQuaternion(cos(ϕ/2), 0, sin(ϕ/2), 0) # rotation ϕ about y
    # rotate particle 1 by q1 only (stays in yz plane)
    # rotate particle 2 by q2, then q3 but in original frame so order swapped
    rotations = [q1, q3*q2]
    Cluster(positions, rotations, sizes, [material for _∈ 1:2], [type for _∈ 1:2])
end



"""
    cluster_helix(N, a, b, c, R, Λ, δ = π/4, δ_0 = 0, handedness="left",
        material = "Au", type="particle")

Helical cluster of N identical particles with axis along z
- `N`: number of particles
- `a,b,c`: semi-axes along x,y,z
- `R`: helix radius
- `Λ`: helix pitch
- `δ`: angle between subsequent particles
- `δ`: starting angle
- `handedness`: "left" or "right"
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_helix(4, 20, 20, 30, 50, 300)
```

"""
function cluster_helix(N, a, b, c, R, Λ, δ = π/4, δ_0 = 0, handedness="left",
    material = "Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles

    s = handedness == "left" ? -1 : +1
    ϕ = collect(1:N) * δ .+ δ_0
    x = R * cos.(ϕ)
    y = R * sin.(ϕ)
    z = s * ϕ * Λ/(2π)
    z .-= s*maximum(s * z)/2 # if <0, add, else remove

    # euler angles calculation
    x′ =  - y
    y′ =   x
    z′ =  s * Λ / (2π)
    n = @. sqrt(x′^2 + y′^2 + z′^2)

    φ =  atan.(y′, x′)
    θ =  acos.(z′ ./ n)
    ψ = 0.0 # don't care for axisymmetric particles

    positions = SVector.(x,y,z)
    # rotations = UnitQuaternion.(RotZYZ.(φ, θ, ψ))
    rotations = UnitQuaternion.(inv.(RotZYZ.(φ, θ, ψ)))

    Cluster(positions, rotations, sizes, [material for _∈ 1:N], [type for _∈ 1:N])
end



"""
    cluster_line(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

Line of N identical particles in the x direction
- `N`: number of particles (approximate if not exact square)
- `a,b,c`: semi-axes along x,y,z
- `Λ`: array pitch
- `φ, θ, ψ`: particle Euler angles
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_line(3, 500, 20, 20, 30, 0, 0, 0)
```

"""
function  cluster_line(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    rotations = [UnitQuaternion(RotZYZ(φ, θ, ψ)) for ii in 1:N] # identical particles

    positions = SVector.(-(N-1)*Λ/2:Λ:(N-1)*Λ/2, 0.0, 0.0)

    Cluster(positions, rotations, sizes, [material for _∈ 1:N], [type for _∈ 1:N])
end



"""
    cluster_array(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

Square array of N identical particles in the xy plane
- `N`: number of particles (approximate if not exact square)
- `a,b,c`: semi-axes along x,y,z
- `Λ`: array pitch
- `φ, θ, ψ`: particle Euler angles
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_array(4, 500, 20, 20, 30, 0, 0, 0)
```

"""
function  cluster_array(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

    N′ = floor(sqrt(N))
    N = N′^2 # actual number,  may have fewer than original N particles

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    rotations = [UnitQuaternion(RotZYZ(φ, θ, ψ)) for ii in 1:N] # identical particles

    x =  -(N′-1)*Λ/2:Λ:(N′-1)*Λ/2
    positions = SVector.(Iterators.product(x, x, 0.0))[:]

    Cluster(positions, rotations, sizes, [material for _∈ 1:N], [type for _∈ 1:N])
end
