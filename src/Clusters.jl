
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

    "types::Vector{String}"
    types::Vector{String}

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
    Cluster(positions, rotations, sizes, [material], [type])
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
    Cluster(positions, rotations, sizes, [material for _ ∈ 1:2], [type for _ ∈ 1:2])
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
    Cluster(positions, rotations, sizes, [material for _ ∈ 1:4], [type for _ ∈ 1:4])
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

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], [type for _ ∈ 1:N])
end



"""
    cluster_chain(N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle")

Line of N identical particles in the x direction
- `N`: number of particles (approximate if not exact square)
- `a,b,c`: semi-axes along x,y,z
- `Λ`: array pitch
- `φ, θ, ψ`: particle Euler angles
- `material`: String referencing the material of every particle
- `type`: String, "point" dipole or "particle"

# Examples

```
cluster_chain(3, 500, 20, 20, 30, 0, 0, 0)
```

"""
function cluster_chain(N, Λ, a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    rotations = [inv(QuatRotation(RotZYZ(α, β, γ))) for _ ∈ 1:N] # identical particles

    positions = SVector.(-(N - 1)*Λ/2:Λ:(N-1)*Λ/2, zero(eltype(Λ)), zero(eltype(Λ)))

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], [type for _ ∈ 1:N])
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

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], [type for _ ∈ 1:N])
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


function sample_random(N)

    α = π * (2 * rand(N) .- 1) # uniform [-pi,pi]
    β = acos.(2 * rand(N) .- 1) # cos-uniform [-1,1]

    positions = @. SVector(cos(α) * sin(β), sin(α) * sin(β), cos(β))
    return positions

end

# hard core pseudo-random
function sample_random_hc(N, exclusion; maxiter=1e3, k=30)


    if sqrt(4 * pi * 1^2 / N) < (pi * exclusion^2)
        @warn "The requested number of points will not fit"
    end

    #initial sample
    s = sample_random(N + k)

    indices = trues(N + k) # all assumed good

    for ii in 1:(N+k) # points to test
        for jj in (ii+1):(N+k)
            dist = norm(s[ii] - s[ii])
            if (dist < exclusion)  # this ii point is bad
                indices[ii] = false
                break # bad point, no need to test further
            end
        end
    end

    todo = (sum(indices) < N)
    # if >= than N, we're done, return N positions 
    if !todo
        @info "no iteration needed"
        pick = findall(indices)
        s = s[pick[1:N]]
        return s
    end

    # otherwise, replace bad and try again
    iter = 0
    while todo
        bad = findall(.!indices)
        p = sum(.!indices)
        @info "iteration $iter, $p bad points so far"
        s[bad] = sample_random(p) # replace bad ones with new random

        for i in eachindex(bad) # points to test
            badi = bad[i]
            indices[badi] = true # assume it is good until shown otherwise
            for j in 1:(N+k)
                if badi == j
                    continue
                end
                dist = norm(s[badi] - s[j])
                if dist < exclusion  # this new point is bad
                    indices[badi] = false
                    break # no need to test further
                end
            end
        end
        # if more than N, we're done
        if sum(indices) >= N
            @info "iterations successful"
            # info "success, %i points out of %i generated", sum(indices), N)
            pick = findall(indices)
            s = s[pick[1:N]]
            return s
        end
        if iter >= maxiter
            @warn "max number of iterations reached"
        end

        todo = (sum(indices) < N) && (iter < maxiter)
        iter = iter + 1
    end

    # we should not get here ideally
    return s

end

# spherical interpolation
# https://en.wikipedia.org/wiki/Slerp
function slerp(p0, p1, d)
    R = norm(p0)
    dt = d / R # angle corresponding to desired spacing
    p0n = p0 / R
    p1n = p1 / norm(p1)
    Omega = acos(dot(p0n, p1n))
    t = dt / Omega
    p = sin((1 - t) * Omega) / sin(Omega) * p0 + sin(t * Omega) / sin(Omega) * p1
    return p
end


function sample_landings(N, threshold_d, dimer_d; maxiter=1e3, k=30)

    if sqrt(4 * pi * 1^2 / N) < (pi * threshold_d^2)
        @warn "The requested number of points will not fit"
    end


    #initial sample
    s = sample_random(N + k)
    sold = s

    indices = trues(N + k) # all assumed good
    dimers = .!indices

    # first pass, checking distances
    for i in 1:(N+k) # points to test
        for j in (i+1):(N+k)
            dist = norm(s[i] - s[j])
            if (dist < threshold_d)  # this i point is bad
                indices[i] = false
                break # bad point, no need to test further
            end
        end
    end

    todo = (sum(indices) < N)

    # if more than N, we're done without dimers, return first N positions 
    if !todo
        @info "no iteration needed, zero dimers"
        return (s=s[1:N], dimers=dimers[1:N])
    end

    # otherwise, move pairs too close to set dimer distance and try again
    iter = 0
    while todo

        number_bad = sum(.!indices)

        @info "iteration $iter, $number_bad bad"
        for i in 1:(N+k) # points to test
            for j in (i+1):(N+k)
                dist = norm(s[i] - s[j])
                if (dist < threshold_d)  # this pair of points is too close, make it a dimer
                    dimers[i] = true
                    dimers[j] = true
                    # now assume this pair is fine until proven otherwise below
                    indices[i] = true
                    indices[i] = true
                    # shift j along the great circle to fixed separation d
                    newp = slerp(s[i], s[j], dimer_d)
                    s[j] = newp
                end
            end
        end

        # redo pass, checking distances
        for i in 1:(N+k) # points to test
            for j in (i+1):(N+k)
                dist = norm(s[i] - s[j])
                if (dist < threshold_d)  # this i point is bad
                    indices[i] = false
                    break # bad point, no need to test further
                end
            end
        end

        number_bad = sum(.!indices)
        @info "now $number_bad bad"
        # if more than N, we're done
        if sum(indices) >= N
            @info "iterations successful"
            pick = findall(indices)
            return (s=s[pick[1:N]], dimers=dimers[pick[1:N]])
        end

        if iter >= maxiter
            @warn "max number of iterations reached"
        end

        todo = (sum(indices) < N) && (iter < maxiter)
        iter = iter + 1
    end

    # we should not get here ideally
    return (s=s[1:N], dimers=dimers[1:N])

end


function cluster_shell_landings(N, a, R, threshold_d=0.5, dimer_d=0.8; monomer_mat="NileBlueM", dimer_mat="NileBlueD", type="point")

    # logic: simulate random landings
    # then check for pairs too close
    # move them about by slerp until no longer chocking
    # while loop until no more collisions
    # tag dimers as such

    s, dimers = sample_landings(N, threshold_d, dimer_d)

    positions = R .* s

    N = length(positions)
    sizes = [SVector(a, a, a) for _ ∈ 1:N] # identical spheres

    # particles isotropic, orientation irrelevant
    trace = a + a + a # whatever scalings we inputed
    sizes = [SVector(trace / 3, trace / 3, trace / 3) for _ ∈ 1:N] # identical particles
    rotations = [@SVector rand(3) for _ ∈ 1:N]

    materials = [monomer_mat for _ ∈ 1:N]
    materials[dimers] .= dimer_mat

    quaternions = [inv(QuatRotation(RotZYZ(r[1], r[2], r[3]))) for r in rotations]

    Cluster(positions, quaternions, sizes, materials, [type for _ ∈ 1:N])
end


function cluster_shell(N, a, b, c, R; orientation="radial", position="fibonacci", material="Rhodamine", type="point", min_exclusion=1.0)

    if position == "fibonacci"
        positions = R .* sample_fibonacci(N)
    elseif position == "random"
        positions = R .* sample_random(N)
    elseif position == "pseudo-random"
        full_area = 4π * R^2
        area_pp = full_area / N
        max_radius = sqrt(area_pp / π)
        max_exclusion = 2 * max_radius  # gives some room
        # exclusion = max_exclusion
        @info "max_exclusion: $max_exclusion"
        exclusion = min(min_exclusion, max_exclusion)
        positions = R .* sample_random_hc(N, exclusion / R)

    else
        @info "No position means fibonacci"
        positions = R .* sample_fibonacci(N)
    end

    N = length(positions) # might be +1
    sizes = [SVector(a, b, c) for _ ∈ 1:N] # identical particles

    if orientation == "radial"

        rotations = map(x -> SVector(atan(x[2], x[1]), acos(x[3] / R), 0), positions)

    elseif orientation == "iso" # make particles isotropic, orientation irrelevant
        trace = a + b + c # whatever scalings we inputed
        sizes = [SVector(trace / 3, trace / 3, trace / 3) for _ ∈ 1:N] # identical particles
        rotations = [@SVector rand(3) for _ ∈ 1:N]

    elseif orientation == "flatiso" # make particles plane-isotropic, orientation radial
        trace = a + b + c # whatever scalings we inputed
        sizes = [SVector(trace / 2, trace / 2, 0) for _ ∈ 1:N] # identical particles
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

    Cluster(positions, quaternions, sizes, [material for _ ∈ 1:N], [type for _ ∈ 1:N])

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

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], [type for _ ∈ 1:N])

end
