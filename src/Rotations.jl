
# splatted version of Rotations.RotZYZ
RotZYZ(v::SVector) = RotZYZ(v...)

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
function euler_active(φ, θ, ψ=0)

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
function euler_passive(φ, θ, ψ=0)

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

euler_unitvector(v::SVector) = euler_unitvector(v[1], v[2]) # ψ irrelevant here



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

    cosθ = cos(θ)
    sinθ = sin(θ)

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

