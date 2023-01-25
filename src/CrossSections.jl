
## cross-sections

struct CrossSections{T}
    extinction::Array{T}
    absorption::Array{T}
    scattering::Array{T}
end

"""
    extinction(kn::Real, P::Array{Complex}, Ein::Array{Complex})

Extinction cross-section for each incident angle

- `kn`: wavenumber in incident medium
- `P`:   `3N_dip x N_inc` matrix, polarisations for all incidences
- `Ein`: `3N_dip x N_inc` matrix, incident field for all incidences

"""
function extinction!(Cext, kn, P, Ein)

    N_inc = size(P, 2)

    for jj in 1:N_inc
        Cext[jj] = 4π * kn * imag(dot(Ein[:, jj], P[:, jj])) # E^* P
    end

    return Cext
end


"""
    absorption(kn::Real, P::Array{Complex}, E::Array{Complex})

Absorption cross-section for each incident angle

- `kn`: wavenumber in incident medium
- `P`: `3N_dip x N_inc` matrix, polarisations for all incidences
- `E`: `3N_dip x N_inc` matrix, total field for all incidences

"""
function absorption!(Cabs, kn, P, E)

    N_inc = size(P, 2)

    for jj in 1:N_inc
        Cabs[jj] = 4π * kn * (imag(dot(E[:, jj], P[:, jj])) -
                              kn^3 * 2 / 3 * real(dot(P[:, jj], P[:, jj])))
    end

    return Cabs
end


# NOTE this seems actually redundant, as boild down to current formula
#
# this is the second version of CDE, where alpha is static
# and radiation reaction is included in the coupling instead
# we convert the results using Eqs 52 and 22 of Markel 2019
# function absorption2!(Cabs, kn, E, P, AlphaBlocks)

#     N_inc = size(P, 2)
#     N_dip = Int(size(P, 1) / 3)
#     Gᵢᵢ = 2im / 3 * k^3 * I # self-reaction
#     # χ: static polarisabilities, we actually need 
#     # local field e = χ⁻¹ p 
#     # e = χ⁻¹α E = (I + Gα) E
#     e = similar(E)
#     # for i in eachindex(AlphaBlocks) # alpha not needed
#     for i in 1:N_dip
#         ii = 3i-2:3i
#         # e = (I + Gᵢᵢ * AlphaBlocks[i]) * E[ii, :]

#         e[ii, :] = E[ii, :] + Gᵢᵢ * P[ii, :] # local field
#         # TODO compute per dipole absorption
#     end

#     for j in 1:N_inc       
#         # Cabs[jj] = 4π * kn * (imag(dot(E[:, jj], P[:, jj])) -
#         # kn^3 * 2 / 3 * real(dot(P[:, jj], P[:, jj])))
#         Cabs[j] = 4π * kn * imag(dot(e[:, j], P[:, j]))
#     end

#     # Cabs = imag(sum(e .* conj(P), dims=1))
#     return Cabs
# end

"""
    scattering(positions::Vector{SVector{3}}, ScatteringVectors::Vector{SVector{3}}, weights::Vector{Real}, kn::Real, P::Array{Complex})

Scattering cross-section for each incident angle, obtained by numerical cubature
over the full solid angle of scattering directions

- `positions`: vector of cluster particle positions
- `ScatteringProjectorz`: `N_inc`-vector of far-field directions
- `weights`: `N_inc`-vector of cubature weights
- `kn`: wavenumber in incident medium
- `P`: `3N_dip x N_inc` matrix, polarisations for all incidences

"""
function scattering!(Csca, positions, ScatteringVectors, weights, kn, P)

    N_dip = length(positions)
    N_inc = size(P, 2)
    N_sca = length(ScatteringVectors)

    # note: maybe this will be needed, though eltype seems to do the trick
    # https://stackoverflow.com/questions/41843949/julia-lang-check-element-type-of-arbitrarily-nested-array
    T = eltype(Csca)
    Isca = zeros(T, N_sca, N_inc) # temp. storage of FF intensities for all scattering directions

    for i = 1:N_sca # loop over scattering angles

        # # unit vector in the scattering direction
        n = ScatteringVectors[i] # rotation of Oz is the third column of Rm

        # far-field "propagator" [kind of]
        nn = n * transpose(n)
        G = (I - nn)

        # temporary storage of net far-field [sum_j Esca[dipole j]]
        # for a given scattering direction
        Esca = zeros(Complex{T}, 3, N_inc)
        for j = 1:N_dip
            jj = (j-1)*3+1:j*3 # find current dipole

            rj = positions[j]
            nrj = dot(n, rj)
            phase = exp(-1im * kn * nrj)

            Esca = Esca + phase * G * P[jj, :]

        end

        # Esca is now the net FF in direction ii
        # Isca[ii, :] = real(sum(Esca .* conj(Esca), dims=1)) # |Esca|^2
        Isca[i, :] = sum(abs2.(Esca), dims=1)
    end

    # now integrate Isca over all scattering angles
    # (for each incident angle)

    Csca[:] = 4π * kn^4 * transpose(weights) * Isca
    return Csca
end




"""
oa_c_ext(kn, G0, B, diagA, shift)

 Computes the OA absorption cross-section from the solution of CDA equations

 PARAMETERS:
 - kn [scalar] wavenumber in medium
 - G 3NrxN_inc matrix
 - W 3NrxN_inc matrix
 - invV 3NrxN_inc matrix

 RETURNS: absorption cross-section analytically orientation-averaged

 DEPENDS: Asym

 FAMILY: low_level, cross_section

"""
function oa_c_ext(kn, G0, B, diagA, shift)
    # TODO experimental, cleanup

    # cext = 2*pi*kn * trace(Asym(G0)*ctranspose(B)*Asym(diagA));
    cext = 2π / kn^2 * tr(shift * Asym(G0) * B' * Asym(diagA'))
    return cext
end


"""
oa_c_abs(k, G0, invA, diagInvAlpha)

 Computes the OA absorption cross-section from the solution of CDA equations

 PARAMETERS:
 - kn [scalar] wavenumber in medium
 - G0 3Nrx3Nr matrix of free-space Green tensor at positions r_i
 - invA 3Nrx3Nr inverse of the interaction matrix
 - diagInvAlpha 3Nrx3Nr inverse of the block-diagonal matrix of polarisabilities

 RETURNS: absorption cross-section analytically orientation-averaged

 DEPENDS: Asym

 FAMILY: low_level, cross_section

"""
function oa_c_abs(k, G0, invA, diagInvAlpha)
    # TODO experimental, cleanup

    # cext = 2*pi*kn * trace(Asym(G0)*ctranspose(B)*Asym(diagA));
    cabs = 2π / k^2 * tr(invA * Asym(G0) * invA' * Asym(diagInvAlpha))
    return cabs
end

# anti-Hermitian part of a matrix
Asym(A) = (A - A') / 2im
