
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
```

"""
function cluster_single(a, b, c, α = 0.0, β = 0.0, γ = 0.0, material = "Au", type="particle")
    sizes = [SVector{3}(a, b, c)]
    positions = [SVector{3}(0.0, 0.0, 0.0)]
    angles = [SVector{3}(α, β, γ)]
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
julia> cluster_dimer(10,1,2,3)
```

"""
function cluster_dimer(d, a, b, c, dihedral = 0.0, material = "Au", type="particle")
    sizes = [SVector{3}(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector{3}(0.0, y, 0.0) for y in (-d/2, d/2)]
    angles = [SVector{3}(0.0, dihedral, 0.0),
              SVector{3}(0.0, 0.0, 0.0)]
    Cluster(positions, angles, sizes, material, type)
end
