# __precompile__()

module CoupledDipole

using LinearAlgebra
using StaticArrays
using Rotations
using SpecialFunctions: besselj, besselh # for Mie
using FastGaussQuadrature: gausslegendre
using DataInterpolations # for tabulated epsilon and alpha
using QuadGK # for ellipsoids, should use GL as well to reduce deps
# convenience
using Makie
using DataFrames
using ProgressMeter

import Rotations: RotZYZ

include("Rotations.jl")
include("Utils.jl")
include("Clusters.jl")
include("Mie.jl")
include("Materials.jl")
include("CrossSections.jl")
include("IncidentField.jl")
include("NearField.jl")
include("HighLevel.jl")
include("PostProcessing.jl")
include("Visual.jl")


# CoupledDipole
export interaction_matrix_labframe!
export interaction_matrix_labframe
export oa_analytical
export polarisation!
export polarisation
export iterate_field!
# IncidentField
export incident_field_pw!
export incident_field_pw
export incident_field_probe
export incident_field_dip!
# Clusters
export Cluster
export cluster_single
export cluster_dimer
export cluster_quadrimer
export cluster_helix
export cluster_chain
export cluster_array
export cluster_shell
export cluster_shell_landings
# CrossSections
export extinction!
export absorption!
export scattering!
export oa_c_ext
export oa_c_abs
# Materials
export Material
export alpha_wrapper
export epsilon_Ag
export epsilon_Au
export epsilon_Si
export epsilon_water
export lorentzian
export alpha_lorentzmolecule
export alpha_rh6g
export alpha_embed
export alpha_scale
export alpha_kuwata
export alpha_particles
export depolarisation_spheroid
export depolarisation_ellipsoid
# Mie
export mie_susceptibility
export ricatti_bessel
export mie_ff
# HighLevel
export spectrum_dispersion
export spectrum_oa
export spectrum_oa_analytical
export map_nf
export scattering_pattern
# PostProcessing
export dispersion_df
export oa_df
export incidence_labels
# Utils
export expand_grid
export pmap_df
export quadrature_lgwt
export cubature_sphere
export euler_active
export euler_passive
export euler_unitvector
export axis_angle
export RotZYZ
export spheroid_ar
export is_inside
export ellipsoid
# Visual
export visualise_makie
export visualise_threejs
# ScatteredField
export local_field
export scattered_field
export far_field

## core functions

"""
    interaction_matrix_labframe!(A,
        kn, R::Vector{SVector{3}},
        AlphaBlocks::Vector{SMatrix{3,3}})

Interaction matrix

- `F`: `3N_dip x 3N_dip` interaction matrix
- `k`: wavenumber in incident medium
- `R`: `N_dip`-vector of 3-Svectors of particle positions
- `AlphaBlocks`: `N_dip`-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function interaction_matrix_labframe!(F, k, R, AlphaBlocks)

    N_dip = length(R)

    # nested loop over dipole pairs
    for i = 1:N_dip
        ii = 3i-2:3i

        for j = i+1:N_dip
            jj = 3j-2:3j

            rⱼ_to_rᵢ = R[i] - R[j]
            rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
            n = rⱼ_to_rᵢ / rᵢⱼ
            nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n
            nx = SMatrix{3,3}(0, n[3], n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p

            expikror = exp(im * k * rᵢⱼ) / rᵢⱼ

            # # EE and MM coupling tensor
            # Aᵢⱼ =
            #     expikror * (
            #         k^2 * (nxn - I) +
            #         (im * k * rᵢⱼ - 1) / (rᵢⱼ^2) * (3 * nxn - I)
            #     )
            kr = k * rᵢⱼ
            k2expikror = k^2 * exp(im * kr) / rᵢⱼ

            Aᵢⱼ = k2expikror * (
                (nxn - I) + (im / kr - 1 / kr^2) * (3 * nxn - I)
            )

            αᵢ = AlphaBlocks[i]
            αⱼ = AlphaBlocks[j]
            Aᵗᵢⱼ = transpose(Aᵢⱼ)

            # assign blocks
            @views F[ii, jj] = Aᵢⱼ * αⱼ
            @views F[jj, ii] = Aᵗᵢⱼ * αᵢ

            # use views to avoid copies (?)
            # Av1 = @view A[ind_jj, ind_kk]
            # Av2 = @view A[ind_kk, ind_jj]
            # fill!(Av1, Ajk * AlphaBlocks[kk])
            # fill!(Av2, transpose(Ajk) * AlphaBlocks[jj])

        end
    end

    return F
end



function interaction_matrix_labframe(k, R, AlphaBlocks)
    N_dip = length(R)
    # what is the type of arrays to initialise?
    proto_r = R[1][1] # position type
    proto_α = AlphaBlocks[1][1, 1] # complex polarisability
    T = typeof(1im * k * proto_r + proto_α) # blocks are ~ exp(ikr) or R * α

    F = Matrix{T}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    interaction_matrix_labframe!(F, k, R, AlphaBlocks)
    return F
end

# anti-Hermitian part of a square matrix
asym(A) = (A - adjoint(A)) / 2im

"""
oa_analytical(kn, cl.positions, AlphaBlocks)

Interaction matrix

- `k`: wavenumber in incident medium
- `R`: `N_dip`-vector of 3-Svectors of particle positions
- `AlphaBlocks`: `N_dip`-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function oa_analytical(k, R, AlphaBlocks)

    N_dip = length(R)
    # what is the type of arrays to initialise?
    proto_r = R[1][1] # position type
    proto_α = AlphaBlocks[1][1, 1] # complex polarisability
    T = typeof(1im * k * proto_r + proto_α) # blocks are ~ exp(ikr) or R * α

    F = Matrix{T}(I, 3N_dip, 3N_dip) # type inferred from cl.positions

    G₀ = 2im / 3 * k^3 * F
    diagAlpha = zeros(T, (3N_dip, 3N_dip))
    diagInvAlpha0 = zeros(T, (3N_dip, 3N_dip))

    # nested loop over dipole pairs
    for i = 1:N_dip
        ii = 3i-2:3i

        αᵢ = AlphaBlocks[i]

        @views diagAlpha[ii, ii] = αᵢ
        # TODO note this won't work for uniaxial alphas
        @views diagInvAlpha0[ii, ii] = inv(αᵢ) + 2im / 3 * k^3 * I

        for j = i+1:N_dip
            jj = 3j-2:3j

            rⱼ_to_rᵢ = R[i] - R[j]
            rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
            n = rⱼ_to_rᵢ / rᵢⱼ
            nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n
            nx = SMatrix{3,3}(0, n[3], n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p

            expikror = exp(im * k * rᵢⱼ) / rᵢⱼ

            kr = k * rᵢⱼ
            k2expikror = k^2 * exp(im * kr) / rᵢⱼ

            Aᵢⱼ = k2expikror * (
                (nxn - I) + (im / kr - 1 / kr^2) * (3 * nxn - I)
            )

            αⱼ = AlphaBlocks[j]
            Aᵗᵢⱼ = transpose(Aᵢⱼ)

            # assign blocks        
            @views G₀[ii, jj] = -Aᵢⱼ
            @views G₀[jj, ii] = -Aᵗᵢⱼ
            @views F[ii, jj] = Aᵢⱼ * αⱼ
            @views F[jj, ii] = Aᵗᵢⱼ * αᵢ

        end
    end

    B = diagAlpha / F
    cabs = -2π / k^2 * real(tr(B * asym(G₀) * adjoint(B) * asym(diagInvAlpha0)))
    cext = 2π / k^2 * real(tr(asym(G₀) * asym(B)))
    return cabs, cext
end


"""
    polarisation!(P, E, AlphaBlocks)

Self-consistent dipole moments from the electric field, P = αE

- `P`: `3N_dip x N_inc` matrix, polarisations for all incidences
- `E`: `3N_dip x N_inc` matrix, total field for all incidences
- `AlphaBlocks`: `N_dip`-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function polarisation!(P, E, AlphaBlocks)

    # loop over N_dip particles
    for i in eachindex(AlphaBlocks)
        ii = 3i-2:3i

        @views P[ii, :] = AlphaBlocks[i] * E[ii, :] # all incidence angles at once
    end

    return P
end

function polarisation(E, AlphaBlocks)

    P = similar(E)
    polarisation!(P, E, AlphaBlocks)

    return P
end


"""
    iterate_field!(E, P, σ_ext, Ein, G, kn, AlphaBlocks, tol=1e-8, maxiter=1000)

Order-of-scattering iteration of the total field

- `E`: `3N_dip x N_inc` matrix, total field for all incidences
- `P`: `3N_dip x N_inc` matrix, polarisations for all incidences
- `Ein`: `3N_dip x N_inc` matrix, total field for all incidences
- `AlphaBlocks`: `N_dip`-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function iterate_field!(
    E,
    P,
    σ_ext,
    Ein,
    F,
    kn,
    AlphaBlocks;
    tol=1e-4,
    maxiter=1000
)

    error = Inf
    niter = 1

    # initialise, first interaction
    Etmp = deepcopy(Ein)
    E[:] = Ein
    polarisation!(P, Ein, AlphaBlocks)
    extinction!(σ_ext, kn, P, Ein)
    σ_prev = deepcopy(σ_ext)

    while (error > tol) && (niter < maxiter)

        Etmp = (I - F) * Etmp
        E = E + Etmp # add new order contribution
        polarisation!(P, E, AlphaBlocks)
        extinction!(σ_ext, kn, P, Ein)
        error = maximum(abs.((σ_ext - σ_prev)) ./ abs.(σ_ext + σ_prev))
        # println(error)
        # println(niter)
        niter += 1


    end

end


end
