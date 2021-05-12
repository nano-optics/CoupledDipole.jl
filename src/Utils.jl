
## utils

"""
    euler_active(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- φ: Euler angle (longitude, in [0,2π])
- θ: Euler angle (colatitude, in [0,2π])
- ψ: Euler angle (rotation around z", in [0,2π])


@example
   euler_active(π/2,0,0)
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

euler_active(v::SVector) = euler_passive(v...)

"""
    euler_passive(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- φ: Euler angle (longitude, in [0,2π])
- θ: Euler angle (colatitude, in [0,2π])
- ψ: Euler angle (rotation around z", in [0,2π])


@example
   euler_passive(π/2,0,0)
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

euler_passive(v::SVector) = euler_passive(v...)

## orientation-averaging

@doc raw"""
    quadrature_lgwt(N::Int, a::Real, b::Real)

N-point Gauss-Legendre quadrature over [a,b] interval
- N: number of nodes
- a,b: bounds

$$
\int_a^b f(x)\,dx=\frac{b-a}{2}\int_{-1}^{1}f\left(\frac{b-a}{2} x + \frac{a+b}{2}\right)\,dx.
$$

@example
   quadrature_lgwt(6,0,3)
"""
function quadrature_lgwt(N, a, b)
    n, w = gausslegendre(N)
    (nodes = (b - a) / 2 * n .+ (a + b) / 2,
     weights = (b - a) / 2 * w)
end

"""
    cubature_sphere(N::Int, method::String)

N-point Gauss-Legendre quadrature over [a,b] interval
- N: number of nodes
- method: cubature method (only 'gl' currently implemented)

Returns a Cubature object containing 2 arrays (Nx3 nodes and Nx1 weights)
Note: using array instead of tuple for weights because we'll probably use them in a scalar product in orientation-averaging
For nodes there is less of a reason, but it can be convenient to visualise the nodes.

The cubature is normalised by 4π such that a unit integrand approximates 1.

@example
   cubature_sphere(36*18)
"""
function cubature_sphere(N, method = "gl")
    #might have slightly more than N total points
    rndN = Integer(ceil(sqrt(N / 2.0)))

    alpha = quadrature_lgwt(2 * rndN, 0, 2π)
    beta = quadrature_lgwt(rndN, 0, 1)

    nodes = [
        SVector{3}(a, acos(2b - 1), 0.0) for
        b in beta.nodes, a in alpha.nodes
    ]
    weights = hcat([a * b for b in beta.weights, a in alpha.weights]...)

    (nodes = nodes[:], weights = 1 / (2π) * weights[:])
end
