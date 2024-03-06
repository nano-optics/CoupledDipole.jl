using StaticArrays
using SpecialFunctions

function ricatti_bessel(x, lmax)
    lnu = 0.5 .+ (0:lmax)
    bj = hcat([besselj(ν, x) * sqrt(pi / 2 * x) for ν in lnu])
    bh = hcat([besselh(ν, x) * sqrt(pi / 2 * x) for ν in lnu])
    ψ′ = [bj[j] - j * bj[j+1] / x for j in 1:lmax]
    ξ′ = [bh[j] - j * bh[j+1] / x for j in 1:lmax]
    ψ = bj[2:lmax+1]
    ξ = bh[2:lmax+1]
    return ψ, ξ, ψ′, ξ′
end


function t_mat(radius, n, μ, λ, l)

    k = 2π * n / λ
    kₛ = 2π * n[1:end-1] ./ λ  # wavenumber j-th
    kₚ = 2π * n[2:end] ./ λ    # wavenumber j+1-th
    xₛ = kₛ .* radius                 # size parameter x = k_j   * r_j
    xₚ = kₚ .* radius                 # size parameter x = k_j+1 * r_j   
    η = n[1:end-1] ./ n[2:end]        # normalized refractive index n_j / n_j+1
    μ = μ[1:end-1] ./ μ[2:end]     # normalized permeability μ_j / μ_j+1
    lmax = length(l)

    # -------------------------------------------------------------------------
    # initialise empty list of list of transfer matrices (one for each interface and each l)
    _z = zeros(SMatrix{2,2})
    tmb = [[_z for i ∈ 1:lmax] for j ∈ eachindex(xₛ)]
    tmf = similar(tmb)
    teb = similar(tmb)
    tef = similar(tmb)

    # -------------------------------------------------------------------------
    ## TRANSFER MATRICES for each interface
    # -------------------------------------------------------------------------
    for j in eachindex(xₛ)

        # compute RB functions for current interface
        # note: should prestore them and extract index j and j+1
        ψₛ, ξₛ, ψₛ′, ξₛ′ = ricatti_bessel(xₛ[j], lmax)
        ψₚ, ξₚ, ψₚ′, ξₚ′ = ricatti_bessel(xₚ[j], lmax)

        for i in 1:lnm
            # note: filled columnwise
            # 1  3
            # 2  4
            tmb[i][j] = -1im * SMatrix{2,2}(
                ξₛ′[i] * ψₚ[i] * η[j] - ξₛ[i] * ψₚ′[i] * μ[j],
                -ψₛ′[i] * ψₚ[i] * η[j] + ψₛ[i] * ψₚ′[i] * μ[j],
                ξₛ′[i] * ξₚ[i] * η[j] - ξₛ[i] * ξₚ′[i] * μ[j],
                -ψₛ′[i] * ξₚ[i] * η[j] + ψₛ[i] * ξₚ′[i] * μ[j])
            # -------------------------------------------------------------------------
            tmf[i][j] = -1im * SMatrix{2,2}(
                ξₚ′[i] * ψₛ[i] / η[j] - ξₚ[i] * ψₛ′[i] / μ[j],
                -ψₚ′[i] * ψₛ[i] / η[j] + ψₚ[i] * ψₛ′[i] / μ[j],
                ξₚ′[i] * ξₛ[i] / η[j] - ξₚ[i] * ξₛ′[i] / μ[j],
                -ψₚ′[i] * ξₛ[i] / η[j] + ψₚ[i] * ξₛ′[i] / μ[j])
            # -------------------------------------------------------------------------         
            teb[i][j] = -1im * SMatrix{2,2}(
                ξₛ′[i] * ψₚ[i] * μ[j] - ξₛ[i] * ψₚ′[i] * η[j],
                -ψₛ′[i] * Sxₚ[i, j] * μ[j] + ψₛ[i] * ψₚ′[i] * η[j],
                ξₛ′[i] * ξₚ[i] * μ[j] - ξₛ[i] * ξₚ′[i] * η[j],
                -ψₛ′[i] * ξₚ[i] * μ[j] + ψₛ[i] * ξₚ′[i] * η[j])
            # -------------------------------------------------------------------------
            tef[i][j] = -1im * SMatrix{2,2}(
                ξₚ′[i] * ψₛ[i] / μ[j] - ξₚ[i] * ψₛ′[i] / η[j],
                -ψₚ′[i] * ψₛ[i] / μ[j] + ψₚ[i] * ψₛ′[i] / η[j],
                ξₚ′[i] * ξₛ[i] / μ[j] - ξₚ[i] * ξₛ′[i] / η[j],
                -ψₚ′[i] * ξₛ[i] / μ[j] + ψₚ[i] * ξₛ′[i] / η[j])
            # -------------------------------------------------------------------------
        end
    end

end