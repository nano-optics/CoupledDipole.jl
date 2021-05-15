
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
- `P`:   3N_dip x N_inc matrix, polarisations for all incidences
- `Ein`: 3N_dip x N_inc matrix, incident field for all incidences

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

- `kn`: wavenumber in incident medium
- `P`: 3N_dip x N_inc matrix, polarisations for all incidences
- `E`: 3N_dip x N_inc matrix, total field for all incidences

"""
function absorption!(Cabs, kn, P, E)

    N_inc = size(P, 2)

    for jj in 1:N_inc
        Cabs[jj] = 4π*kn*(imag(dot(E[:,jj], P[:,jj])) -
                         kn^3 * 2/3 * real(dot(P[:,jj], P[:,jj])))
    end

    return  Cabs
end


"""
    scattering(positions, angles, weights, kn::Real, P::Array{Complex})

Scattering cross-section for each incident angle, obtained by numerical cubature
over the full solid angle of scattering directions

- `positions`: vector of cluster particle positions
- `ScatteringProjectorz`: `N_inc`-vector of 3x3 projector Smatrices
- `weights`: N_inc-vector of cubature weights
- `kn`: wavenumber in incident medium
- `P`: 3N_dip x N_inc matrix, polarisations for all incidences

"""
function scattering!(Csca, positions, ScatteringVectors, weights, kn, P)

    N_dip = length(positions)
    N_inc = size(P, 2)
    N_sca = length(ScatteringVectors)

    # note: maybe this will be needed, though eltype seems to do the trick
    # https://stackoverflow.com/questions/41843949/julia-lang-check-element-type-of-arbitrarily-nested-array
    T = eltype(Csca)
    Isca = zeros(T, N_sca, N_inc); # temp. storage of FF intensities for all scattering directions

    for ii = 1:N_sca # loop over scattering angles

        # # unit vector in the scattering direction
        n = ScatteringVectors[ii] # rotation of Oz is the third column of Rm

        # far-field "propagator" [kind of]
        nn = n*transpose(n)
        G = (I - nn)

        # temporary storage of net far-field [sum_j Esca[dipole j]]
        # for a given scattering direction
        Esca = zeros(Complex{T},3,N_inc)
        for jj = 1:N_dip

            rj = positions[jj]
            nrj = dot(n,rj)
            indjj = (jj-1)*3+1:jj*3 # find current dipole
            phase = exp(-1im*kn*nrj)

            Esca = Esca +  phase * G * P[indjj, :]

        end

        # Esca is now the net FF in direction ii
        Isca[ii,:] = real(sum(Esca.*conj(Esca),dims=1)) # |Esca|^2

    end

    # now integrate Isca over all scattering angles
    # (for each incident angle)

    Csca[:] = 4π * kn^4 * transpose(weights) * Isca
    return Csca
end
