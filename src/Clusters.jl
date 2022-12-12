
## cluster definitions


"""
    Cluster(positions, angles, sizes, material)

Particle cluster specification.

"""
struct Cluster{T1,T2,T3}

    "positions::Vector{SVector{3,T1}}"
    positions::Vector{SVector{3,T1}}

    "rotations::Vector{QuatRotation{T2}}"
    rotations::Vector{QuatRotation{T2}}

    "sizes::Vector{SVector{3,T3}}"
    sizes::Vector{SVector{3,T3}}

    "materials::Vector{String}"
    materials::Vector{String}

    "type::String"
    type::String

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
function cluster_single(a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")
    sizes = [SVector(a, b, c)]
    positions = [SVector(0.0, 0.0, 0.0)]
    # input parameters are Euler angles
    # inverse as passive rotation needed
    rotations = [inv(QuatRotation(Rotations.RotZYZ(α, β, γ)))]
    Cluster(positions, rotations, sizes, [material], type)
end


"""
    cluster_dimer(d::T, a::T, b::T, c::T, ϕ::T = 0.0, α_1::T = 0.0, α_2::T = 0.0) where T <: Real

Particle cluster consisting of 2 identical particles separated along y
- `a,b,c`: semi-axes along x,y,z
- `ϕ`: dihedral angle between both particles seen along the y-axis
- `α_1,α_2`: angle of each particle with the y axis
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_dimer(80, 10, 10, 20)
```

"""
function cluster_dimer(d, a, b, c, ϕ=0.0, α_1=0.0, α_2=0.0, material="Au", type="particle")
    sizes = [SVector(a, b, c) for _ ∈ 1:2] # identical particles
    positions = [SVector(zero(eltype(d)), y, zero(eltype(d))) for y in (-d / 2, d / 2)]
    q1 = QuatRotation(cos(α_1 / 2), sin(α_1 / 2), zero(eltype(α_1)), zero(eltype(α_1))) # rotation α_1 about x
    q2 = QuatRotation(cos(α_2 / 2), sin(α_2 / 2), zero(eltype(α_1)), zero(eltype(α_1))) # rotation α_2 about x
    q3 = QuatRotation(cos(ϕ / 2), zero(eltype(α_1)), sin(ϕ / 2), zero(eltype(α_1))) # rotation ϕ about y
    # rotate particle 1 by q1 only (stays in yz plane)
    # rotate particle 2 by q2, then q3 but in original frame so order swapped
    rotations = inv.([q1, q3 * q2]) # inverse as passive rotation needed
    Cluster(positions, rotations, sizes, [material for _ ∈ 1:2], type)
end


"""
   cluster_quadrimer(r, d, material="Au", type="particle")

Particle cluster consisting of 4 identical spheres in a chiral geometry
- `r`: sphere radius
- `d`: gap between spheres
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_quadrimer(80, 10)
```

"""
function cluster_quadrimer(r, d, material="Au", type="particle")
    sizes = [SVector(r, r, r) for _ ∈ 1:4] # identical particles
    p1 = SVector(0.0, 0.0, 0.0)
    p2 = SVector(2r + d, 0.0, 0.0)
    p3 = SVector(0.0, 2r + d, 0.0)
    p4 = SVector(0.0, 2r + d, 2r + d)
    positions = [p1, p2, p3, p4]
    q = QuatRotation(1.0, 0.0, 0.0, 0.0)
    rotations = [q, q, q, q] # useless since spheres
    Cluster(positions, rotations, sizes, [material for _ ∈ 1:4], type)
end

# function [cl] = (r, s, d, theta)

#     cl.angles = zeros(3, 4)
#     cl.sizes = r * ones(3, 4)
#     p1 = [0; 0; 0]
#     p2 = [s; 0; 0]
#     p3 = [0; 0; d]
#     p4 = [s * cos(theta); s * sin(theta); d]

#     cl.positions = [p1 p2 p3 p4]

# end



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
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_helix(4, 20, 20, 30, 50, 300)
```

"""
function cluster_helix(N, a, b, c, R, Λ, δ, δ_0=0, handedness="left",
    material="Au", type="particle")

    sizes = [SVector(a, b, c) for _ ∈ 1:N] # identical particles

    s = handedness == "left" ? -1 : +1
    ϕ = collect(1:N) * δ .+ δ_0
    x = R * cos.(ϕ)
    y = R * sin.(ϕ)
    z = s * ϕ * Λ / (2π)
    z .-= s * maximum(s * z) / 2 # if < 0, add, else remove

    # euler angles calculation
    x′ = -y
    y′ = x
    z′ = s * Λ / (2π)
    n = @. sqrt(x′^2 + y′^2 + z′^2)

    φ = atan.(y′, x′)
    θ = acos.(z′ ./ n)
    ψ = 0.0 # don't care for axisymmetric particles

    positions = SVector.(x, y, z)
    # rotations = QuatRotation.(RotZYZ.(φ, θ, ψ))
    # inverse as passive rotation needed
    rotations = QuatRotation.(inv.(RotZYZ.(φ, θ, ψ)))

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)
end



"""
    cluster_line(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

Line of N identical particles in the x direction
- `N`: number of particles (approximate if not exact square)
- `a,b,c`: semi-axes along x,y,z
- `Λ`: array pitch
- `φ, θ, ψ`: particle Euler angles
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_line(3, 500, 20, 20, 30, 0, 0, 0)
```

"""
function cluster_line(N, Λ, a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    rotations = [inv(QuatRotation(RotZYZ(α, β, γ))) for _ ∈ 1:N] # identical particles

    positions = SVector.(-(N - 1)*Λ/2:Λ:(N-1)*Λ/2, zero(eltype(Λ)), zero(eltype(Λ)))

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)
end



"""
    cluster_array(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

Square array of N identical particles in the xy plane
- `N`: number of particles (approximate if not exact square)
- `a,b,c`: semi-axes along x,y,z
- `Λ`: array pitch
- `φ, θ, ψ`: particle Euler angles
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_array(4, 500, 20, 20, 30, 0, 0, 0)
```

"""
function cluster_array(N, Λ, a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")

    N′ = floor(sqrt(N))
    N = N′^2 # actual number,  may have fewer than original N particles

    sizes = [SVector(a, b, c) for _ ∈ 1:N] # identical particles
    rotations = [inv(QuatRotation(RotZYZ(α, β, γ))) for _ ∈ 1:N] # identical particles

    x = -(N′ - 1)*Λ/2:Λ:(N′-1)*Λ/2
    positions = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)
end

function sample_fibonacci(N)

    N = ifelse(mod(N, 2) == 0, N + 1, N)

    j = 0:N-1
    β = @. acos(1 - (2 * j + 1) / N)
    ϕ = (1 + √5) / 2
    α = @. (2π * j / ϕ) % 2π

    positions = @. SVector(cos(α) * sin(β), sin(α) * sin(β), cos(β))
    return positions

end

function cluster_shell(N, a, b, c, R; orientation="radial", material="Rhodamine", type="point")

    positions = R .* sample_fibonacci(N)
    N = length(positions) # might be +1
    sizes = [SVector(a, b, c) for _ ∈ 1:N] # identical particles

    if orientation == "radial"

        rotations = map(x -> SVector(atan(x[2], x[1]), acos(x[3] / R), 0), positions)

    elseif orientation == "flat"

        # strategy a bit suboptimal but doesn't matter
        # first find the normal vector
        # then create two basis vectors in the tangent plane
        # sample a random linear combination of the two as our flat-random direction

        φ = map(x -> atan(x[2], x[1]), positions)
        θ = map(x -> acos(x[3] / R), positions)
        τ_1 = map(φ -> SVector(-sin(φ), cos(φ), 0), φ)
        τ_2 = map((φ, θ) -> SVector(cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)), φ, θ)
        τ = map((τ_1, τ_2) -> rand(1) .* τ_1 .+ rand(1) .* τ_2, τ_1, τ_2)
        rotations = map(τ -> SVector(atan(τ[2], τ[1]),
                acos(τ[3] / sqrt(sum(τ .^ 2))), 0), τ)

    else
        @info "No orientation means random orientation"
        rotations = [@SVector rand(3) for _ ∈ 1:N]

    end

    quaternions = [inv(QuatRotation(RotZYZ(r[1], r[2], r[3]))) for r in rotations]

    Cluster(positions, quaternions, sizes, [material for _ ∈ 1:N], type)

end



function cluster_ball(N, a, R; material="AuCM", type="point")
    #example:  N=4169 -> exact for R=10
    # possible improvement would be to use the inscribed cube
    # but also only do 1/8 of the solid angle by symmetry, etc.
    # cube inscribed
    # s = 2a/√3

    # Gauss approx formula suggests
    # Ng = π^(3 / 2) * 4 / (3√π) * R^3
    # so the radius that gives us ≈ N points is 
    Rg = Int(ceil((N / (π^(3 / 2) * 4 / (3√π)))^(1 / 3)))
    Rg2 = Rg^2
    s = R / Rg
    # integer lattice points within that approximate ball, rescaled to actual ball
    positions = [SVector(s * x, s * y, s * z) for x in -N:N for y in -N:N for z in -N:N if x^2 + y^2 + z^2 <= Rg2]
    N = length(positions) # should be similar to N requested
    @info "$N points generated for that ball"
    sizes = [SVector(a, a, a) for _ ∈ 1:N] # identical particles
    rotations = [inv(QuatRotation(1, 0, 0, 0.0)) for _ ∈ 1:N] # no rotations, spheres anyway...

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)

end
