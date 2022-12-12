"""
scattered_field(λ, probes, positions, dipoles)

- `λ`: wavelength
- `probes`: SVector of spatial positions where the field is to be evaluated
- `positions`: SVector of dipole positions
- `dipoles`: 3Ndipx2Ninc matrix of self-consistent dipole moments

"""
function scattered_field(k, Rdip, Rpro, P)

    N_dip = length(Rdip)
    N_pro = length(Rpro)
    N_inc = size(P, 2)
    Z₀ = 376.730313668 # free-space impedance

    # Esca = zeros(eltype(P), (3N_pro, N_inc))
    # actually easier to work with array of vectors
    Esca = [@SVector(zeros(eltype(P), (3))) for _ ∈ 1:N_inc*N_pro]
    Bsca = [@SVector(zeros(eltype(P), (3))) for _ ∈ 1:N_inc*N_pro]

    # for loop over N_probe locations
    for i = 1:N_pro
        # store contribution at current probe location for all directions
        Etmp = zeros(eltype(P), (3, N_inc))
        Btmp = zeros(eltype(P), (3, N_inc))

        for j = 1:N_dip
            jj = 3j-2:3j
            # source to probe
            rⱼ_to_rᵢ = Rpro[i] - Rdip[j]
            rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
            n = rⱼ_to_rᵢ / rᵢⱼ
            nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n
            nx = SMatrix{3,3}(0, n[3], n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p

            expikror = exp(im * k * rᵢⱼ) / rᵢⱼ

            # EE 
            Aᵢⱼ =
                expikror * (
                    k^2 * (nxn - I) +
                    (im * k * rᵢⱼ - 1) / (rᵢⱼ^2) * (3 * nxn - I)
                )

            # EM 
            Bᵢⱼ = expikror * (nx - I) * (k^2 + im * k / rᵢⱼ)

            # contribution from j-th source dipole, for all incidences
            Etmp += Aᵢⱼ * P[jj, :]
            Btmp += Bᵢⱼ * P[jj, :]

        end
        # store scattered field 3-vector at probe location for each incidence
        for l ∈ 1:N_inc
            Esca[N_inc*(i-1)+l] = Etmp[:, l]
            Bsca[N_inc*(i-1)+l] = Z₀ .* Btmp[:, l]
        end

    end

    return Esca
end

