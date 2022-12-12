"""
scattered_field(λ, probes, positions, dipoles)

- `λ`: wavelength
- `probes`: SVector of spatial positions where the field is to be evaluated
- `positions`: SVector of dipole positions
- `dipoles`: 3Ndipx2Ninc matrix of self-consistent dipole moments

"""
function scattered_field(kn, Rdip, Rpro, P)

    N_dip = length(Rdip)
    N_pro = length(Rpro)
    N_inc = size(P, 2)

    # Esca = zeros(eltype(P), (3N_pro, N_inc))
    # actually easier to work with array of vectors
    Esca = [@SVector(zeros(eltype(P), (3))) for _ ∈ 1:N_inc*N_pro]
    Bsca = [@SVector(zeros(eltype(P), (3))) for _ ∈ 1:N_inc*N_pro]

    # for loop over N_probe locations
    for jj = 1:N_pro
        # when Esca was a matrix
        # ind_jj = 3jj-2:3jj
        # where to store the fields in Esca
        ind2_jj = N_inc*(jj-1)+1:N_inc*jj
        Etmp = zeros(eltype(P), (3, N_inc)) # store contribution at current probe location for all directions
        Btmp = zeros(eltype(P), (3, N_inc))
        # for loop over N_dip dipoles
        for kk = 1:N_dip
            ind_kk = 3kk-2:3kk

            # source to probe
            rk_to_rj = Rpro[jj] - Rdip[kk]
            rjk = norm(rk_to_rj, 2)
            rjkhat = rk_to_rj / rjk
            rjkrjk = rjkhat * transpose(rjkhat)

            Ge_jk =
                exp(im * kn * rjk) / rjk * (
                    kn * kn * (rjkrjk - I) +
                    (im * kn * rjk - 1.0) / (rjk * rjk) * (3 * rjkrjk - I)
                )
            Gm_jk =
                exp(im * kn * rjk) / rjk * (
                    kn * kn * (rjkrjk - I) +
                    (im * kn * rjk - 1.0) / (rjk * rjk) * (3 * rjkrjk - I)
                )
            # Esca[ind_jj, :] += Gjk * P[ind_kk, :]
            # now a list
            Etmp += Ge_jk * P[ind_kk, :]
            Btmp += Gm_jk * P[ind_kk, :]

        end
        # store scattered electric field 3-vector at probe location for each incidence
        for i ∈ 1:N_inc
            Esca[N_inc*(jj-1)+i] = Etmp[:, i]
            Bsca[N_inc*(jj-1)+i] = Btmp[:, i]
        end
        # Esca[ind2_jj] = [Etmp[:, l] for l in 1:N_inc]
        # Bsca[ind2_jj] = [Btmp[:, l] for l in 1:N_inc]

    end

    return Esca
end

