
## cross-sections

struct CrossSections{T}
    extinction::Array{T}
    absorption::Array{T}
    scattering::Array{T}
end

"""
    extinction(kn::Real, P::Array{Complex}, Ein::Array{Complex})

Extinction cross-section for each incident angle

- kn: wavenumber in incident medium
- P:   3N_dip x N_inc matrix, polarisations for all incidences
- Ein: 3N_dip x N_inc matrix, incident field for all incidences

"""
function extinction!(Cext, kn, P, Ein)

    N_inc = size(P,2)

    for jj in 1:N_inc
        Cext[jj] = 4π*kn*imag(dot(Ein[:,jj], P[:,jj])) # E^* P
    end

    return  Cext
end


"""
    absorption(kn::Real, P::Array{Complex}, E::Array{Complex})

Absorption cross-section for each incident angle

- kn: wavenumber in incident medium
- P: 3N_dip x N_inc matrix, polarisations for all incidences
- E: 3N_dip x N_inc matrix, total field for all incidences

"""
function absorption!(Cabs, kn, P, E)

    N_inc = size(P, 2)

    for jj in 1:N_inc
        Cabs[jj] = 4π*kn*imag(dot(E[:,jj], P[:,jj]) -
        kn^3 * 2/3 * real(dot(P[:,jj],P[:,jj])))
    end

    return  Cabs
end


# note: need also Csca as numerical cubature, for consistency check with Ext-Abs

function scattering!(Csca, positions, angles, weights, kn, P)

    N_dip = length(positions)
    N_inc = size(P, 2)
    N_sca = length(angles)
    # G = zeros(3, 3)
    # I3 = I
    #
    T = typeof(Csca)
    Isca = zeros(N_sca, N_inc); # temp. storage of FF intensities for all scattering directions

    for ii = 1:N_sca # loop over scattering angles

        # # unit vector in the scattering direction
        Rm = rotation_euler_active(angles[ii]...)
        n = Rm[:,3] # rotation of Oz is the third column of Rm

        # # FF propagator [akin to]
        nn = n*transpose(n)
        G = (I - nn)

        # temporary storage of net far-field [sum_j Esca[dipole j]]
        # for a given scattering direction
        # Esca = zeros(3, N_inc);
        for jj = 1:N_dip

            rj = positions[:,jj]
            nrj = dot(n,rj)
            indjj = (jj-1)*3+1:jj*3 # find current dipole
            phase = exp(-1i*kn*nrj)
            Esca = Esca +  phase * G * P[indjj, :]

        end

        # Esca is now the net FF in direction ii
        Isca[ii,:] = real(sum(Esca.*conj(Esca))) # |Esca|^2

    end

    # need to integrate Isca over all scattering angles
    # (for each incident angle)

    Csca = 4π * kn^4 * transpose(weights) * Isca
    return Csca
end
