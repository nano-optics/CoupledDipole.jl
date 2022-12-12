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

    Esca = zeros(eltype(P), (3N_pro, N_inc))

    # for loop over N_probe locations
    for jj = 1:N_pro
        ind_jj = 3jj-2:3jj

        # for loop over N_dip dipoles
        for kk = 1:N_dip
            ind_kk = 3kk-2:3kk

            # source to probe
            rk_to_rj = Rpro[jj] - Rdip[kk]
            rjk = norm(rk_to_rj, 2)
            rjkhat = rk_to_rj / rjk
            rjkrjk = rjkhat * transpose(rjkhat)

            Gjk =
                exp(im * kn * rjk) / rjk * (
                    kn * kn * (rjkrjk - I) +
                    (im * kn * rjk - 1.0) / (rjk * rjk) * (3 * rjkrjk - I)
                )

            Esca[ind_jj, :] += Gjk * P[ind_kk, :]

        end
    end

    return Esca
end

