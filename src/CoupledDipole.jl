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
export propagator_freespace_labframe!
export incident_field!
export polarisation!
export iterate_field!
# Clusters
export Cluster
export cluster_single
export cluster_dimer
export cluster_helix
export cluster_line
export cluster_array
export cluster_shell
# CrossSections
export extinction!
export absorption!
export scattering!
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
export alpha_blocks!
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
    propagator_freespace_labframe!(A,
        kn, R::Vector{SVector{3}},
        AlphaBlocks::Vector{SMatrix{3,3}})

Interaction matrix

- `A`: `3N_dip x 3N_dip` interaction matrix
- `kn`: wavenumber in incident medium
- `R`: `N_dip`-vector of 3-Svectors of particle positions
- `AlphaBlocks`: `N_dip`-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function propagator_freespace_labframe!(A, kn, R, AlphaBlocks)

    N_dip = length(R)

    # Note: might want to try BlockArray.jl
    # Construct a BlockArray from blocks
    # mortar((A,B),(C,D))
    # would simply indexing, at the expense of another dep

    # nested for loop over N_dip dipoles
    for jj = 1:N_dip
        ind_jj = 3jj-2:3jj

        for kk = (jj+1):N_dip
            ind_kk = 3kk-2:3kk

            rk_to_rj = R[jj] - R[kk]
            rjk = norm(rk_to_rj, 2)
            rjkhat = rk_to_rj / rjk
            rjkrjk = rjkhat * transpose(rjkhat)

            Ajk =
                exp(im * kn * rjk) / rjk * (
                    kn * kn * (rjkrjk - I) +
                    (im * kn * rjk - 1.0) / (rjk * rjk) * (3 * rjkrjk - I)
                )

            # assign blocks
            A[ind_jj, ind_kk] = Ajk * AlphaBlocks[kk]
            A[ind_kk, ind_jj] = transpose(Ajk) * AlphaBlocks[jj]
            # use views to avoid copies (?)
            # Av1 = @view A[ind_jj, ind_kk]
            # Av2 = @view A[ind_kk, ind_jj]
            # fill!(Av1, Ajk * AlphaBlocks[kk])
            # fill!(Av2, transpose(Ajk) * AlphaBlocks[jj])

        end
    end

    return A
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
