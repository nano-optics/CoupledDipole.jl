
## material functions

struct Material{T}
    wavelength::Vector{T}
    medium::Dict{String,Function}
end


"""
    epsilon_Ag(λ::T) where T <: Real

Drude model for the dielectric function of silver in the visible region
- λ: wavelength in nm

@example
   epsilon_Ag(632.8)
"""
function epsilon_Ag(λ)
    4 * (1.0 - 1.0 / (282.0^2 * (1 / λ^2 + im / (17000 * λ))))
end

"""
    epsilon_Au(λ::T) where T <: Real

Drude model for the dielectric function of gold in the visible region
- λ: wavelength in nm

@example
   epsilon_Au(632.8)
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
        exp(-1im * φ) / (1 / λ1 + 1 / λ + 1im / μ_p2)
    )
end


"""
    lorentzian(λ::T, α_k::T, λ_k::T, µ_k::T) where T <: Real

Complex Lorentz function, to describe complex polarisabilities
- λ: wavelength in nm
- α_k: oscillator strength in S.I. units
- λ_k: oscillator wavelength in nm
- µ_k: damping in S.I. units

@example
   lorentzian(632.8)
"""
function lorentzian(λ, α_k = 5.76e-38, λ_k = 526.0, µ_k = 1.0e4)
    -α_k * λ_k / µ_k *
    (1.0 - 1.0 / (1.0 - (λ_k / λ)^2 - 1im * (λ_k^2 / (µ_k * λ))))
end

"""
    alpha_bare(λ::T, α_∞::T, α_k::Array{T}, λ_k::Array{T}, µ_k::Array{T}) where T <: Real

Complex scalar polarisability, as sum of lorentz oscillators
- λ: wavelength in nm
- α_k: oscillator strength(s) in S.I. units
- λ_k: oscillator wavelength(s) in nm
- µ_k: damping(s) in S.I. units

Default values mimic the main resonance of Rhodamine 700

@example
   alpha_bare(632.8)
"""
function alpha_bare(λ, α_∞ = 9.6e-39, α_k = 5.76e-38, λ_k = 526.0, µ_k = 1.0e4)

    ε_0 = 8.8541878128e-12
    nm3 = 1e27

    α = α_∞
    for kk = 1:length(α_k)
        α += lorentzian(λ, α_k[kk], λ_k[kk], µ_k[kk])
    end
    prefact = nm3 / (4π * ε_0)
    prefact * α
end


"""
    alpha_embedded(α::Complex{T}, medium::T) where T <: Real

Complex scalar polarisability, as sum of lorentz oscillators
- α: bare polarisabilty
- medium: refractive index of embedding medium

Default values mimic the main resonance of Rhodamine 700

@example
   alpha_embedded(alpha_bare(632.8))
"""
function alpha_embedded(α, medium = 1.33)
    ε_m = medium^2
    L = (ε_m + 2) / 3
    1 / ε_m * L^2 * α
end


"""
    alpha_rescale_molecule(alpha, sizes::Vector{SVector{3}})

Principal polarisability components of a particle, rescaled along each principal axis
- α: scalar polarisabilty
- sizes: array of 3-vectors to scale along each principal axis

@example
   sizes = [SVector{3}(1.0, 2.0, 3.0) for i in 1:4]
   alpha_rescale_molecule(alpha_bare(632.8), sizes)
"""
function alpha_rescale_molecule(alpha, sizes)
    alpha .* (sizes ./ sum.(sizes))
end


"""
    depolarisation_spheroid(a, b, c)

Depolarisation factor of a spheroid
- a semi-axis along x and y
- b semi-axis along x and y (unused)
- c semi-axis along z

@example
   depolarisation_spheroid(1, 1, 1.5)
"""
function depolarisation_spheroid(a, b, c)
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


"""
    alpha_kuwata(wavelength, epsilon, V, axis, L, medium)

Principal polarisability components of a spheroidal particle
- λ wavelength
- ε complex dielectric function
- Size SVector with 3 semi-axes of the spheroid
- χ depolarisation factor
- ε_m dielectric function of surrounding medium

@example
   alpha_kuwata(500, -10+1im, SVector(30, 30, 50), 1.33^2)
"""
function alpha_kuwata(λ, ε, Size, ε_m)

    V = 4 / 3 * π * prod(Size)
    x_0 = @. 2 * π * Size / λ
    χ = depolarisation_spheroid(Size...)

    A = @. -0.4865 * χ - 1.046 * χ^2 + 0.8481 * χ^3
    B = @. 0.01909 * χ + 0.19999 * χ^2 + 0.6077 * χ^3

    denom = @. (χ + ε_m / (ε - ε_m)) + A * ε_m * x_0^2 + B * ε_m^2 * x_0^4 -
       1im / 3 * 4 * π^2 * ε_m^(3 / 2) * V / λ^3

    return ((V / (4 * π)) ./ denom)
end



"""
    alpha_spheroids(λ, ε, ε_m, Sizes)

Principal polarisability components of N spheroidal particles
- λ wavelength
- ε complex dielectric function
- ε_m dielectric function of surrounding medium
- Sizes: Vector of SVectors of particle sizes

@example
   alpha_spheroids(500, -10+1im, 1.33^3, [SVector(30, 30, 50) for i in 1:4])
"""
function alpha_spheroids(λ, ε, ε_m, Sizes)
    return (map(s -> alpha_kuwata(λ, ε, s, ε_m), Sizes))
end
