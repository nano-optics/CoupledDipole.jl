
## utils

"""
    euler_active(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- `φ`: Euler angle (longitude, in [0,2π])
- `θ`: Euler angle (colatitude, in [0,2π])
- `ψ`: Euler angle (rotation around z", in [0,2π])

# Examples

```jldoctest
julia> round.(euler_active(π/2,0,0), digits=5)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.0  -1.0  0.0
  1.0   0.0  0.0
 -0.0   0.0  1.0
```

"""
function euler_active(φ, θ, ψ = 0)

    cosφ = cos(φ)
    cosψ = cos(ψ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    sinψ = sin(ψ)
    sinθ = sin(θ)

    # filled by columns
    R = SMatrix{3,3}(
        cosφ * cosθ * cosψ - sinφ * sinψ,
        sinφ * cosθ * cosψ + cosφ * sinψ,
        -sinθ * cosψ,
        -cosφ * cosθ * sinψ - sinφ * cosψ,
        -sinφ * cosθ * sinψ + cosφ * cosψ,
        sinθ * sinψ,
        cosφ * sinθ,
        sinφ * sinθ,
        cosθ,
    )

    return R

end

# splatted version
euler_active(v::SVector) = euler_active(v...)

"""
    euler_passive(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- `φ`: Euler angle (longitude, in [0,2π])
- `θ`: Euler angle (colatitude, in [0,2π])
- `ψ`: Euler angle (rotation around z", in [0,2π])

# Examples

```jldoctest
julia> round.(euler_passive(π/2,0,0), digits=5)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.0  1.0  -0.0
 -1.0  0.0   0.0
  0.0  0.0   1.0
```

"""
function euler_passive(φ, θ, ψ = 0)

    cosφ = cos(φ)
    cosψ = cos(ψ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    sinψ = sin(ψ)
    sinθ = sin(θ)

    # filled by columns
    R = SMatrix{3,3}(
        cosφ * cosθ * cosψ - sinφ * sinψ,
        -cosφ * cosθ * sinψ - sinφ * cosψ,
        cosφ * sinθ,
        sinφ * cosθ * cosψ + cosφ * sinψ,
        -sinφ * cosθ * sinψ + cosφ * cosψ,
        sinφ * sinθ,
        -sinθ * cosψ,
        sinθ * sinψ,
        cosθ,
    )

    return R

end

# splatted version
euler_passive(v::SVector) = euler_passive(v...)




"""
    euler_unitvector(φ::Real, θ::Real)

Unit vector along direction φ, θ

- `φ`: Euler angle (longitude, in [0,2π])
- `θ`: Euler angle (colatitude, in [0,2π])

# Examples

```jldoctest
julia> euler_unitvector(π/2, 0)
3-element SVector{3, Float64} with indices SOneTo(3):
 0.0
 0.0
 1.0
```

"""
function euler_unitvector(φ, θ)

    cosφ = cos(φ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    sinθ = sin(θ)
    n = SVector(cosφ * sinθ, sinφ * sinθ, cosθ)

    return n

end

euler_unitvector(v::SVector) = euler_unitvector(v[1],v[2]) # ψ irrelevant here



@doc raw"""
    axis_angle(v = SVector(0, 1, 0), θ)

3D rotation matrix from axis-angle representation

- `v`: SVector
- `θ`: rotation angle

# Examples

```
julia> axis_angle(SVector(0, 1, 0), π/4)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.707107  0.0  0.707107
  0.0       1.0  0.0
 -0.707107  0.0  0.707107
```

"""
function axis_angle(v, θ)

    cosθ = cos(θ); sinθ = sin(θ)

    R = SMatrix{3,3}(
        cosθ + v[1]^2 * (1 - cosθ),
        v[2] * v[1] * (1 - cosθ) + v[3] * sinθ,
        v[3] * v[1] * (1 - cosθ) - v[2] * sinθ,
        #
        v[1] * v[2] * (1 - cosθ) - v[3] * sinθ,
        cosθ + v[2]^2 * (1 - cosθ),
        v[3] * v[2] * (1 - cosθ) + v[1] * sinθ,
        #
        v[1] * v[3] * (1 - cosθ) + v[2] * sinθ,
        v[2] * v[3] * (1 - cosθ) - v[1] * sinθ,
        cosθ + v[3]^2 * (1 - cosθ)
    )

    return R
end


## orientation-averaging

@doc raw"""
    quadrature_lgwt(N::Int, a::Real, b::Real)

N-point Gauss-Legendre quadrature over [a,b] interval
- `N`: number of nodes
- `a,b`: bounds

```math
\int_a^b f(x)\,dx=\frac{b-a}{2}\int_{-1}^{1}f\left(\frac{b-a}{2} x + \frac{a+b}{2}\right)\,dx.
```

# Examples

```
julia> quadrature_lgwt(6,0,3)
(nodes = [0.10129572869527204, 0.5081859203006033, 1.1420712208752046, 1.8579287791247954, 2.4918140796993966, 2.898704271304728], weights = [0.25698673856875537, 0.5411423595722078, 0.7018709018590369, 0.7018709018590369, 0.5411423595722078, 0.25698673856875537])
```

"""
function quadrature_lgwt(N, a, b)
    n, w = gausslegendre(N)
    (nodes = (b - a) / 2 * n .+ (a + b) / 2,
     weights = (b - a) / 2 * w)
end

"""
    cubature_sphere(N::Int, method::String)

N-point cubature on the sphere
- `N`: number of nodes
- `method`: cubature method (only 'gl' currently implemented)

Returns a Cubature object containing 2 arrays (N'x3 nodes and N'x1 weights), N'≈N
Note: using array instead of tuple for weights because
we'll use them in a scalar product in orientation-averaging
For nodes there is less of a reason, but it can be convenient to visualise the nodes.

The cubature is normalised by 4π such that a unit integrand approximates 1.

# Examples

```
julia> cubature_sphere(6)
(nodes = SVector{3, Float64}[[0.43625314334650644, 2.1862760354652844, 0.0], [0.43625314334650644, 0.9553166181245093, 0.0], [2.0735107047038173, 2.1862760354652844, 0.0], [2.0735107047038173, 0.9553166181245093, 0.0], [4.2096746024757685, 2.1862760354652844, 0.0], [4.2096746024757685, 0.9553166181245093, 0.0], [5.846932163833079, 2.1862760354652844, 0.0], [5.846932163833079, 0.9553166181245093, 0.0]], weights = [0.08696371128436346, 0.08696371128436346, 0.16303628871563655, 0.16303628871563655, 0.16303628871563655, 0.16303628871563655, 0.08696371128436346, 0.08696371128436346])
```

"""
function cubature_sphere(N, method = "gl")
    #might have slightly more than N total points
    rndN = Integer(ceil(sqrt(N / 2.0)))

    φ = quadrature_lgwt(2rndN, 0, 2π) # longitude
    cθ = quadrature_lgwt(rndN, 0, 1)   # cos(colatitude)

    nodes = [
        SVector{3}(a, acos(2b - 1), 0.0) for
        b in cθ.nodes, a in φ.nodes
    ]
    weights = hcat([a * b for b in cθ.weights, a in φ.weights]...)

    (nodes = nodes[:], weights = 1 / (2π) * weights[:])
end
