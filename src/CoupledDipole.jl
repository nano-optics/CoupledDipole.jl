# __precompile__()

module CoupledDipole

using LinearAlgebra
using StaticArrays
using Rotations
using SpecialFunctions: besselj, besselh # for Mie
using FastGaussQuadrature: gausslegendre
using DataInterpolations # for tabulated epsilon and alpha
using QuadGK # for ellipsoids, should use GL as well to reduce deps
using Makie
using DataFrames

import Rotations: RotZYZ

include("Rotations.jl")
include("Utils.jl")
include("Clusters.jl")
include("Mie.jl")
include("Materials.jl")
include("CrossSections.jl")
include("NearField.jl")
include("HighLevel.jl")
include("PostProcessing.jl")
include("Visual.jl")


# CoupledDipole
export interaction_matrix_labframe!
export incident_field!
export polarisation!
export iterate_field!
# Clusters
export Cluster
export cluster_single
export cluster_dimer
export cluster_quadrimer
export cluster_helix
export cluster_chain
export cluster_array
export cluster_shell
# CrossSections
export extinction!
export absorption!
export scattering!
export oa_c_ext
# Materials
export Material
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
# PostProcessing
export dispersion_df
export oa_df
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
# Visual
export visualise_makie
export visualise_threejs
# NearField
export scattered_field

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

            # EE and MM coupling tensor
            Aᵢⱼ =
                expikror * (
                    k^2 * (nxn - I) +
                    (im * k * rᵢⱼ - 1) / (rᵢⱼ^2) * (3 * nxn - I)
                )

            αᵢ = AlphaBlocks[i]
            αⱼ = AlphaBlocks[j]
            Aᵗᵢⱼ = transpose(Aᵢⱼ)

            # assign blocks
            F[ii, jj] = Aᵢⱼ * αⱼ
            F[jj, ii] = Aᵗᵢⱼ * αᵢ

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

"""
    incident_field!(Ein,
        Ejones,
        kn, R::Vector{SVector{3}},
        IncidenceRotations)

Incident field at particle positions

- `Ein`: `3N_dip x N_inc` matrix, right-hand side of coupled-dipole system
- `Ejones`: tuple of 2 2-Svectors defining 2 orthogonal Jones polarisations
- `kn`: wavenumber in incident medium
- `R`: `N_dip`-vector of 3-Svectors of particle positions
- `IncidenceRotations`: `N_inc`-vector of rotation 3-Smatrices

"""
function incident_field!(Ein, Ejones, kn, R, IncidenceRotations)

    N_inc = length(IncidenceRotations)
    N_dip = length(R)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

    for jj in eachindex(IncidenceRotations)
        Rm = IncidenceRotations[jj]
        # Rm * [0;0;1] == 3rd column
        k_hat = kn * Rm[:, 3]
        E1_r = Rm * Evec1
        E2_r = Rm * Evec2
        for kk = 1:N_dip
            kR = dot(k_hat, R[kk]) # everything real

            # use views to avoid copies (?)
            # Ev1 = @view Ein[kk*3-2:kk*3, jj]
            # Ev2 = @view Ein[kk*3-2:kk*3, jj+N_inc]
            # fill!(Ev1, E1_r * exp(im * kR))
            # fill!(Ev2, E2_r * exp(im * kR))

            Ein[kk*3-2:kk*3, jj] = E1_r * exp(im * kR)
            Ein[kk*3-2:kk*3, jj+N_inc] = E2_r * exp(im * kR)
        end
    end
    return Ein
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
    for ii in eachindex(AlphaBlocks)
        ind = 3ii-2:3ii
        # use views to avoid copies (?)
        # Pv = @view P[ind, :]
        # Ev = @view E[ind, :]
        # fill!(Pv, AlphaBlocks[ii] * Ev)

        P[ind, :] = AlphaBlocks[ii] * E[ind, :] # all incidence angles at once
    end

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
