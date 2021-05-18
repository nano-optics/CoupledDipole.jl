
## cluster definitions

# NOTE: for now all particles have the same material
# easy to extend if needed, by making it a vector
# and matching a dictionary in HighLevel functions
struct Cluster{T1,T2,T3}
    positions::Vector{SVector{3,T1}}
    angles::Vector{SVector{3,T2}}
    sizes::Vector{SVector{3,T3}}
    material::String
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

```jldoctest
julia> cluster_single(1.0,2.0,3.0)
Cluster{Float64, Float64, Float64}(SVector{3, Float64}[[0.0, 0.0, 0.0]], SVector{3, Float64}[[0.0, 0.0, 0.0]], SVector{3, Float64}[[1.0, 2.0, 3.0]], "Au", "particle")
```

"""
function cluster_single(a, b, c, α = 0.0, β = 0.0, γ = 0.0, material = "Au", type="particle")
    sizes = [SVector(a, b, c)]
    positions = [SVector(0.0, 0.0, 0.0)]
    angles = [SVector(α, β, γ)]
    Cluster(positions, angles, sizes, material, type)
end


"""
    cluster_dimer(d::T, a::T, b::T, c::T, dihedral::T = 0.0, α_1::T = 0.0, α_2::T = 0.0) where T <: Real

Particle cluster consisting of 2 identical particles separated along y
- `a,b,c`: semi-axes along x,y,z
- `dihedral`: angle between both particles seen along the y-axis
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```jldoctest
julia> cluster_dimer(10, 1, 2, 3)
Cluster{Float64, Float64, Int64}(SVector{3, Float64}[[0.0, -5.0, 0.0], [0.0, 5.0, 0.0]], SVector{3, Float64}[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], SVector{3, Int64}[[1, 2, 3], [1, 2, 3]], "Au", "particle")
```

"""
function cluster_dimer(d, a, b, c, dihedral = 0.0, material = "Au", type="particle")
    sizes = [SVector(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector(0.0, y, 0.0) for y in (-d/2, d/2)]
    angles = [SVector(0.0, dihedral, 0.0),
              SVector(0.0, 0.0, 0.0)]
    Cluster(positions, angles, sizes, material, type)
end



"""
    cluster_helix(N, a, b, c, R, λ, δ = π/4, δ_0 = 0, handedness="left",
        material = "Au", type="particle")

Helical cluster of N identical particles with axis along z
- `N`: number of particles
- `a,b,c`: semi-axes along x,y,z
- `R`: helix radius
- `λ`: helix pitch
- `δ`: angle between subsequent particles
- `δ`: starting angle
- `handedness`: "left" or "right"
- `material`: String referencing the material
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_helix(5, 20, 20, 30, 50, 300)
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

    # angles calculation
    x′ =  - y
    y′ =   x
    z′ =  s * Λ / (2π)
    n = @. sqrt(x′^2 + y′^2 + z′^2)

    φ =  atan.(y′, x′)
    θ =  acos.(z′ ./ n)
    ψ = 0.0 # don't care for axisymmetric particles

    positions = SVector.(x,y,z)
    angles = SVector.(φ,θ,ψ)

    Cluster(positions, angles, sizes, material, type)
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
cluster_line(10, 500, 20, 20, 30, 0, 0, 0)
```

"""
function  cluster_line(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    angles = [SVector(φ, θ, ψ) for ii in 1:N] # identical particles

    positions = SVector.(-(N-1)*Λ/2:Λ:(N-1)*Λ/2), 0.0, 0.0)

    Cluster(positions, angles, sizes, material, type)
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
cluster_array(10, 500, 20, 20, 30, 0, 0, 0)
```

"""
function  cluster_array(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

    N′ = floor(sqrt(N)) # may have fewer than N particles

    sizes = [SVector(a, b, c) for ii in 1:N′^2] # identical particles
    angles = [SVector(φ, θ, ψ) for ii in 1:N′^2] # identical particles

    x =  -(N′-1)*Λ/2:Λ:(N′-1)*Λ/2
    positions = map(x -> SVector(x), Iterators.product(x, x, 0.0))[:]

    Cluster(positions, angles, sizes, material, type)
end
