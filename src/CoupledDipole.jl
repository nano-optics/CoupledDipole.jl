# __precompile__()

module CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature: gausslegendre

#using BlockArrays
#using Base.Threads
#using DataFrames
#using VegaLite
#using Gadfly


include("Clusters.jl")
include("Materials.jl")
include("Utils.jl")
include("CrossSections.jl")
include("HighLevel.jl")

export rotation_matrices!
export alpha_blocks!
export propagator_freespace_labframe!
export incident_field!
export polarisation!
export Cluster
export cluster_single
export cluster_dimer
export extinction!
export absorption!
export scattering!
export Material
export epsilon_Ag
export epsilon_Au
export lorentzian
export alpha_bare
export alpha_kuwata
export alpha_blocks!
export alpha_embedded
export alpha_spheroids
export alpha_rescale_molecule
export depolarisation_spheroid
export quadrature_lgwt
export euler_active
export euler_passive
export cubature_sphere
export spectrum_dispersion
export spectrum_oa

## core functions


# this is not needed, just broadcast
function rotation_matrices!(Rotations, Angles)

    N_dip = length(Angles)

    # loop over N_dip particles
    for ii = 1:N_dip
        Rotations[ii] = euler_passive(Angles[ii]...)
    end

    return Rotations
end


"""
    alpha_blocks(AlphaBlocks::Vector{SMatrix{3,3}},
                 Alpha::Vector{SVector{3}},
                 Angles::Vector{SVector{3}})

Polarisability blocks for the interaction matrix

- AlphaBlocks: N_dip-vector of 3x3 Smatrices (polarisability tensors in the lab frame)
- Alpha: N_dip-vector of principal polarisability components for each particle
- Angles: Euler angles for rotating each particle in the lab frame

"""
function alpha_blocks!(AlphaBlocks, Alpha, Angles)

    N_dip = length(Angles)

    # loop over N_dip particles
    for ii = 1:N_dip
        Rm = euler_passive(Angles[ii]...)
        AlphaBlocks[ii] = Rm' * diagm(Alpha[ii]) * Rm
        # transpose(Rot) * diag(Alpha(ind)) * Rot;
    end

    return AlphaBlocks
end

"""
    propagator_freespace_labframe!(A,
        kn, R::Vector{SVector{3}},
        AlphaBlocks::Vector{SMatrix{3,3}})

Interaction matrix

- A: 3N_dip x 3N_dip interaction matrix
- kn: wavenumber in incident medium
- R: N_dip-vector of 3-Svectors of particle positions
- AlphaBlocks: N_dip-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

# Note: might want/need to try BlockArray.jl
# Construct a BlockArray from blocks
# mortar((A,B),(C,D))

"""
function propagator_freespace_labframe!(A, kn, R, AlphaBlocks)

    N_dip = length(R)

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

        end
    end

    return A
end



# function incident_field!(Ein, Ejones, kn, R, Incidence)
#
#     N_inc = length(Incidence)
#     N_dip = length(R)
#
#     Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
#     Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector
#     Evec_r = similar(Evec1)
#
#     for jj in eachindex(Incidence)
#         Rm = euler_active(Incidence[jj]...)
#         # Rot * [0;0;1] == 3rd column
#         k_hat = kn * transpose(Rm[:, 3])
#         for kk = 1:N_dip
#             Ein[kk*3-2:kk*3, jj] = (Rm * (Evec1 * exp(im * k_hat * R[kk])))
#             Ein[kk*3-2:kk*3, jj+N_inc] = (Rm * (Evec2 * exp(im * k_hat * R[kk])))
#         end
#     end
#     return Ein
# end


"""
    incident_field(Ein,
        Ejones,
        kn, R::Vector{SVector{3}},
        IncidenceRotations)

Incident field at particle positions

- Ein: 3N_dip x N_inc matrix, right-hand side of coupled-dipole system
- Ejones: tuple of 2 2-Svectors defining 2 orthogonal Jones polarisations
- kn: wavenumber in incident medium
- R: N_dip-vector of 3-Svectors of particle positions
- IncidenceRotations: N_inc-vector of 3-Smatrices of rotations

"""
function incident_field!(Ein, Ejones, kn, R, IncidenceRotations)

    N_inc = length(IncidenceRotations)
    N_dip = length(R)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

    for jj in eachindex(IncidenceRotations)
        Rm = IncidenceRotations[jj]
        # Rot * [0;0;1] == 3rd column
        k_hat = kn * Rm[:, 3]
        for kk = 1:N_dip
            kR = dot(k_hat, R[kk])
            Ein[kk*3-2:kk*3, jj] = (Rm * (Evec1 * exp(im * kR)))
            Ein[kk*3-2:kk*3, jj+N_inc] = (Rm * (Evec2 * exp(im * kR)))
        end
    end
    return Ein
end



"""
    polarisation!(P, E, AlphaBlocks)

Incident field at particle positions

- P: 3N_dip x N_inc matrix, polarisations for all incidences
- E: 3N_dip x N_inc matrix, total field for all incidences
- AlphaBlocks: N_dip-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function polarisation!(P, E, AlphaBlocks)

    # loop over N_dip particles
    for ii in eachindex(AlphaBlocks)
        ind = 3ii-2:3ii
        P[ind, :] = AlphaBlocks[ii] * E[ind, :] # all incidence angles at once
    end

    return P
end



end
