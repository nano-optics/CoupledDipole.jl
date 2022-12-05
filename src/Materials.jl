
## material functions

"""
    Material(wavelength, media)

Wavelength-dependent optical properties.

"""
struct Material{T}
    wavelengths::Vector{T}
    media::Dict{String,Function}
end


"""
    epsilon_Ag(λ::Real)

Drude model for the dielectric function of silver in the visible region
- `λ`: wavelength in nm

# Examples

```jldoctest
julia> round(epsilon_Ag(632.8), digits=5)
-16.11377 + 0.74871im
```

"""
function epsilon_Ag(λ)
    4 * (1.0 - 1.0 / (282.0^2 * (1 / λ^2 + 1im / (17000 * λ))))
end

"""
    epsilon_Au(λ::Real)

Extended Drude model for the dielectric function of gold in the visible region
- `λ`: wavelength in nm

# Examples

```jldoctest
julia>  round(epsilon_Au(632.8), digits=5)
-11.40185 + 1.18679im
```

"""
function epsilon_Au(λ)

    ε_∞ = 1.54
    λ_p = 177.5
    μ_p = 14500.0
    A1 = 1.27
    λ1 = 470.0
    μ_p1 = 1900.0
    A2 = 1.1
    λ2 = 325.0
    μ_p2 = 1060.0
    φ = -π / 4

    ε_∞ * (1 - 1 / (λ_p^2 * ((1 / λ)^2 + 1im / (μ_p * λ)))) +
    A1 / λ1 * (
        exp(1im * φ) / (1 / λ1 - 1 / λ - 1im / μ_p1) +
        exp(-1im * φ) / (1 / λ1 + 1 / λ + 1im / μ_p1)
    ) +
    A2 / λ2 * (
        exp(1im * φ) / (1 / λ2 - 1 / λ - 1im / μ_p2) +
        exp(-1im * φ) / (1 / λ2 + 1 / λ + 1im / μ_p2)
    )
end


"""
    lorentzian(λ::Real, α_k::Real, λ_k::Real, µ_k::Real)

Complex Lorentz function, to describe polarisabilities
- `λ`: wavelength in nm
- `α_k`: oscillator strength in S.I. units
- `λ_k`: oscillator wavelength in nm
- `µ_k`: damping in S.I. units

# Examples

```jldoctest
julia> round(lorentzian(632.8)*1e39, digits=5)
6.58095 + 1.35961im
```

"""
function lorentzian(λ, α_k=5.76e-38, λ_k=526.0, µ_k=1.0e4)
    -α_k * λ_k / µ_k *
    (1.0 - 1.0 / (1.0 - (λ_k / λ)^2 - 1im * (λ_k^2 / (µ_k * λ))))
end

"""
    alpha_baremolecule(λ::T, α_∞::T, α_k::Array{T}, λ_k::Array{T}, µ_k::Array{T}) where T <: Real

Complex scalar polarisability, as sum of lorentz oscillators
- `λ`: wavelength in nm
- `α_k`: oscillator strength(s) in S.I. units
- `λ_k`: oscillator wavelength(s) in nm
- `µ_k`: damping(s) in S.I. units

Default values mimic the main resonance of Rhodamine 700

# Examples


```jldoctest
julia> round(alpha_baremolecule(632.8), digits=5)
0.14543 + 0.01222im
```

"""
function alpha_baremolecule(λ, α_∞=9.6e-39, α_k=5.76e-38, λ_k=526.0, µ_k=1.0e4)

    ε₀ = 8.8541878128e-12
    nm3 = 1e27

    α = α_∞
    for kk = eachindex(α_k)
        α += lorentzian(λ, α_k[kk], λ_k[kk], µ_k[kk])
    end
    prefact = nm3 / (4π * ε₀)
    prefact * α
end




"""
    alpha_embed(α::Complex{T}, medium::T) where T <: Real

Effective point polarisability in medium, rescaled by local field correction
- `α`: bare polarisabilty
- `medium`: refractive index of embedding medium

Default values mimic the main resonance of Rhodamine 700

# Examples

```jldoctest
julia> round(alpha_embed(alpha_baremolecule(632.8)), digits=5)
0.12976 + 0.0109im
```

"""
function alpha_embed(α, medium=1.33)
    ε_m = medium^2
    L = (ε_m + 2) / 3
    1 / ε_m * L^2 * α
end


"""
    alpha_scale(alpha, sizes::SVector{3})

Principal polarisability components of a particle, rescaled along each principal axis
- `α`: scalar polarisabilty
- `sizes`: 3-vector to scale along each principal axis

"""
function alpha_scale(alpha, sizes)
    alpha .* (sizes / sum(sizes))
end


"""
alpha_cm(ε, ε_m, Size)

Clausius-Mossotti polarisability
- `ε`: complex dielectric function
- `ε_m`: dielectric function of surrounding medium
- `Size`: particle size, assumed to be spherical
"""
function alpha_cm(ε, ε_m, Size)
    return Size[1]^3 * (ε - ε_m) / (ε + 2 * ε_m)
end


"""
    alpha_particles(λ, ε, ε_m, Sizes)

Principal polarisability components of N particles
- `λ`: wavelength
- `ε`: complex dielectric function
- `ε_m`: dielectric function of surrounding medium
- `Sizes`: Vector of 3-SVectors of particle sizes

# Examples

```
julia> alpha_particles(500, -10+1im, 1.33^3, [SVector(30, 30, 50) for i in 1:4])
4-element Vector{SVector{3, ComplexF64}}:
 [83399.81975161123 + 64172.28157743772im, 83399.81975161123 + 64172.28157743772im, -86034.64340475321 + 79773.72029581512im]
 [83399.81975161123 + 64172.28157743772im, 83399.81975161123 + 64172.28157743772im, -86034.64340475321 + 79773.72029581512im]
 [83399.81975161123 + 64172.28157743772im, 83399.81975161123 + 64172.28157743772im, -86034.64340475321 + 79773.72029581512im]
 [83399.81975161123 + 64172.28157743772im, 83399.81975161123 + 64172.28157743772im, -86034.64340475321 + 79773.72029581512im]
```

"""
function alpha_particles(Epsilon, Sizes, ε_m, λ; prescription="kuwata")

    if prescription == "kuwata"
        return (map((e, s) -> alpha_kuwata(λ, e, ε_m, s), Epsilon, Sizes))
    elseif prescription == "majic"
        return (map((e, s) -> alpha_majic(λ, e, ε_m, s), Epsilon, Sizes))
    elseif prescription == "mie"
        return (map((e, s) -> alpha_mie(λ, e, ε_m, s), Epsilon, Sizes))
    else
        @warning "unknown prescription $prescription"
        return (map((e, s) -> alpha_kuwata(λ, e, ε_m, s), Epsilon, Sizes))
    end
end

"""
ALPHA_MIE exact polarisability from Mie theory

 Given a particle radius, returns the
 wavelength-dependent polarisability

 PARAMETERS:
 - wavelength vector of length Nl
 - epsilon complex vector of length Nl
 - medium refractive index of surrounding medium
 - a radius

 RETURNS: vector of length Nl

 DEPENDS: mini_GDAB [private]

 FAMILY: low_level, polarizability

"""
function alpha_mie(λ, ε, ε_m, Size)

    n_medium = sqrt(ε_m)
    s = sqrt(ε) / n_medium
    k = n_medium * 2π / λ
    x = k * Size[1]
    conversion = 2 / 3 * 1i * k^3
    Γ, Δ, A, B = mie_susceptibility(x, s, 1)

    return Δ / conversion
end

"""
alpha_majic returns a 3NrxNl matrix of principal polarizability components for a spheroid

 Principal polarizability components describing Nr particles of spheroidal shape
 with defined sizes, and using a wavelength-dependent permittivity of length Nl

 PARAMETERS:
 - wavelength
 - epsilon
 - medium
 - sizes

 RETURNS: 3NrxNl matrix of principal polarizability components
 describing spheroids

 FAMILY: user_level, polarizability
"""
function alpha_majic(λ, ε, ε_m, Size)

    εᵣ = (ε / ε_m)
    k = 2π / λ * √ε_m

    a, b, c = Size
    V = 4 / 3 * π * a * b * c

    Lx(e) = -1 / (2 * e^2) * ((1 - e^2) / e * atanh(e) - 1)
    Lz(e) = (1 - e^2) / e^2 * (atanh(e) / e - 1)

    AlphaS(L) = V / 4 / π * (εᵣ - 1) / (1 + (εᵣ - 1) * L)

    if (c ≈ a) & (b ≈ a) # sphere

        L1 = 1 / 3
        e = 0
        Omega1 = 1 / 5 * (εᵣ - 2) / (1 + (εᵣ - 1) * L1)
        AlphaS₁ = AlphaS(L1)
        A1 = AlphaS₁ / (1 - Omega1 * k^2 * c^2 - 1im * 2 / 3 * k^3 * AlphaS₁)
        return (SVector(A1, A1, A1))

    elseif (a ≈ b) # spheroid along z, i.e. standard orientation

        e = sqrt(c^2 - a^2) / c
        L1 = Lx(e)
        L2 = Lz(e)
        Omega1 = 1 / 5 * (εᵣ - 2 + 3 * e^2) / (1 + (εᵣ - 1) * L1) - 12 / 25 * e^2
        Omega2 = 1 / 5 * (εᵣ - 2 - εᵣ * e^2) / (1 + (εᵣ - 1) * L2) + 9 / 25 * e^2
        AlphaS₁ = AlphaS(L1)
        AlphaS₂ = AlphaS(L2)
        A1 = AlphaS₁ / (1 - Omega1 * k^2 * c^2 - 1im * 2 / 3 * k^3 * AlphaS₁)
        A2 = AlphaS₂ / (1 - Omega2 * k^2 * c^2 - 1im * 2 / 3 * k^3 * AlphaS₂)
        return (SVector(A1, A1, A2))

    elseif (b ≈ c) # spheroid along x
        # c -> a 
        # a -> b 
        e = sqrt(a^2 - b^2) / a
        # unchanged below #
        L1 = Lx(e)
        L2 = Lz(e)
        Omega1 = 1 / 5 * (εᵣ - 2 + 3 * e^2) / (1 + (εᵣ - 1) * L1) - 12 / 25 * e^2
        Omega2 = 1 / 5 * (εᵣ - 2 - εᵣ * e^2) / (1 + (εᵣ - 1) * L2) + 9 / 25 * e^2
        AlphaS₁ = AlphaS(L1)
        AlphaS₂ = AlphaS(L2)
        # /unchanged #

        A1 = AlphaS₁ / (1 - Omega1 * k^2 * a^2 - 1im * 2 / 3 * k^3 * AlphaS₁)
        A2 = AlphaS₂ / (1 - Omega2 * k^2 * a^2 - 1im * 2 / 3 * k^3 * AlphaS₂)
        return (SVector(A2, A1, A1)) # swap z and x

    elseif (a ≈ c) # spheroid along y
        # c -> b 
        # a -> a 
        e = sqrt(b^2 - a^2) / b
        # unchanged below #
        L1 = Lx(e)
        L2 = Lz(e)
        Omega1 = 1 / 5 * (εᵣ - 2 + 3 * e^2) / (1 + (εᵣ - 1) * L1) - 12 / 25 * e^2
        Omega2 = 1 / 5 * (εᵣ - 2 - εᵣ * e^2) / (1 + (εᵣ - 1) * L2) + 9 / 25 * e^2
        AlphaS₁ = AlphaS(L1)
        AlphaS₂ = AlphaS(L2)
        # /unchanged #

        A1 = AlphaS₁ / (1 - Omega1 * k^2 * a^2 - 1im * 2 / 3 * k^3 * AlphaS₁)
        A2 = AlphaS₂ / (1 - Omega2 * k^2 * b^2 - 1im * 2 / 3 * k^3 * AlphaS₂)
        return (SVector(A1, A2, A1))# swap z and y

    else # general ellipsoid, not implemented
        @warn "not implemented"
    end

end

"""
    alpha_kuwata(λ, ε, ε_m, Size)

Principal polarisability components of a spheroidal particle
- `λ`: wavelength
- `ε`: complex dielectric function
- `ε_m`: dielectric function of surrounding medium
- `Size`: SVector with 3 semi-axes of the spheroid

# Examples

```
julia> alpha_kuwata(500, -10+1im, 1.33^2, SVector(30, 30, 50))
3-element SVector{3, ComplexF64} with indices SOneTo(3):
  77076.04648078184 + 26235.664281642235im
  77076.04648078184 + 26235.664281642235im
 -98187.15974124733 + 205835.30299929058im
```

"""
function alpha_kuwata(λ, ε, ε_m, Size)

    V = 4π / 3 * prod(Size)
    x = @. 2π / λ * Size
    χ = depolarisation_ellipsoid(Size...)

    A = @. -0.4865 * χ - 1.046 * χ^2 + 0.8481 * χ^3
    B = @. 0.01909 * χ + 0.19999 * χ^2 + 0.6077 * χ^3

    denom = @. (χ + ε_m / (ε - ε_m)) + A * ε_m * x^2 + B * ε_m^2 * x^4 -
               1im / 3 * 4 * π^2 * ε_m^(3 / 2) * V / λ^3

    return ((V / (4π)) ./ denom)
end

"""
  depolarisation_ellipsoid(a, b, c)

Depolarisation factor of an ellipsoid
- `a`: semi-axis along x 
- `b`: semi-axis along y
- `c`: semi-axis along z

# Examples

```jldoctest
julia> round.(depolarisation_ellipsoid(1, 1, 1.5), digits=5)
3-element SVector{3, Float64} with indices SOneTo(3):
 0.38351
 0.38351
 0.23298
```

"""
function depolarisation_ellipsoid(a, b, c)
    if (c ≈ a) & (b ≈ a) # sphere
        Lx = 1 / 3.0
        Ly = Lx
        Lz = Lx
    elseif (a ≈ b) # spheroid along z
        Lx, Ly, Lz = depolarisation_spheroid(a, c)
    elseif (b ≈ c) # spheroid along x
        Lz, Ly, Lx = depolarisation_spheroid(b, a) # swap z and x
    elseif (a ≈ c) # spheroid along y
        Lx, Lz, Ly = depolarisation_spheroid(a, b)  # swap z and y
    else # general case (no closed form)
        V = a * b * c
        Lx = QuadGK.quadgk(q -> 1.0 / ((1.0 + q) * sqrt((q + 1.0) * (q + (b / a)^2) * (q + (c / a)^2))), 0.0, Inf)[1] * V / 2.0 / a^3
        Ly = QuadGK.quadgk(q -> 1.0 / ((1.0 + q) * sqrt((q + 1.0) * (q + (a / b)^2) * (q + (c / b)^2))), 0.0, Inf)[1] * V / 2.0 / b^3
        Lz = QuadGK.quadgk(q -> 1.0 / ((1.0 + q) * sqrt((q + 1.0) * (q + (a / c)^2) * (q + (b / c)^2))), 0.0, Inf)[1] * V / 2.0 / c^3
    end

    return SVector(Lx, Ly, Lz)
end

"""
    depolarisation_spheroid(a, c)

Depolarisation factor of a spheroid
- `a`: semi-axis along x and y
- `c`: semi-axis along z

# Examples

```jldoctest
julia> round.(depolarisation_spheroid(1, 1.5), digits=5)
3-element SVector{3, Float64} with indices SOneTo(3):
 0.38351
 0.38351
 0.23298
```

"""
function depolarisation_spheroid(a, c)
    if (c == a)# sphere
        Lz = 1 / 3
    end

    if (c > a) # prolate
        e = sqrt(1 - (a / c)^2)
        Lz = (1 - e^2) / e^3 * (-e + atanh(e))
    end

    if (c < a) # oblate
        e = sqrt(1 - (c / a)^2)
        Lz = 1 / e^2 * (1 - sqrt(1 - e^2) / e * asin(e))
    end

    Lx = (1 - Lz) / 2

    return (SVector(Lx, Lx, Lz))
end
