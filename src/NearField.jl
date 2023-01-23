
internal_field(p, V, εₘ, εᵣ) = 4π * εₘ / (εᵣ - 1.0) * p / V

"""
scattered_field(probe, k, n, positions, sizes, rotations, P)

- `probe`: SVector, point where the scattered EM field is to be evaluated
- `k`: wavenumber
- `n_medium`: refractive index of medium 
- `positions`: SVector of dipole positions
- `sizes`: SVector of dipole sizes
- `rotations`: SVector of particle rotation matrices
- `P`: 3Ndipx2Ninc matrix of self-consistent dipole moments
- `Epsilon`: list of dielectric function inside each particle (for internal fields)
- `evaluate_inside`: logical, whether to compute the (unphysical) scattered field inside particles

"""
function scattered_field(probe, k, n_medium, positions, sizes, rotations, P, Epsilon=missing; evaluate_inside=true)

    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / 376.730313668 # H = Y₀ E
    c₀ = 299792458 # m/s
    ε₀ = 8.8541878128e-12

    N_dip = length(positions)
    N_inc = size(P, 2)

    Esca = zeros(eltype(P), (3, N_inc))
    Bsca = zeros(eltype(P), (3, N_inc))

    inside, id = is_inside(probe, positions, sizes, rotations)

    if inside & !evaluate_inside # use the dipole moment to approximate Einside
        j = id[]
        jj = 3j-2:3j
        pmod = sum(abs2.(P[jj, :]), dims=1)
        # @info "inside particle $j, $pmod"
        εᵣ = Epsilon[j]
        Esca = internal_field(P[jj, :], 4π / 3 * prod(sizes[j]), n_medium^2, εᵣ)
    end

    if !inside | evaluate_inside # we're outside particles
        for j = 1:N_dip
            jj = 3j-2:3j

            # source to probe
            rⱼ_to_rᵢ = probe - positions[j]
            rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
            n = rⱼ_to_rᵢ / rᵢⱼ
            nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n
            nx = SMatrix{3,3}(0, n[3], -n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p

            expikror = exp(im * k * rᵢⱼ) / rᵢⱼ

            # EE 
            Aᵢⱼ =
                expikror * (
                    k^2 * (I - nxn) -
                    (1 - im * k * rᵢⱼ) / (rᵢⱼ^2) * (I - 3 * nxn)
                )

            # EM 
            Bᵢⱼ = expikror * nx * (k^2 + im * k / rᵢⱼ)

            # contribution from j-th source dipole, for all incidences
            Esca += Aᵢⱼ * P[jj, :]
            Bsca += Bᵢⱼ * P[jj, :]

        end
    end
    # B = Z ε₀ε Gm P̄ 
    # TODO clarify this scaling factor n?
    scale = Z₀ * ε₀ * n_medium
    return Esca, scale * Bsca, inside, id
end

"""
far_field(direction, k, n_medium, cl.positions, cl.sizes,
    ParticleRotations, P, Epsilon)

- `direction`: SVector, point where the scattered EM field is to be evaluated
- `k`: wavenumber
- `n_medium`: refractive index of medium 
- `positions`: SVector of dipole positions
- `P`: 3Ndipx2Ninc matrix of self-consistent dipole moments

"""
function far_field(direction, k, positions, dipoles)

    # Z₀ = 376.730313668 # free-space impedance
    # Y₀ = 1 / 376.730313668 # H = Y₀ E
    # c₀ = 299792458 # m/s
    # ε₀ = 8.8541878128e-12

    N_dip = length(positions)

    Esca = zeros(eltype(dipoles), (3, 3)) # 3 orthogonal sources considered

    # # unit vector in the scattering direction
    n = euler_unitvector(direction)

    # far-field "propagator" [kind of]
    nn = n * transpose(n)
    G = (I - nn)

    # storage of net far-field [sum_j Esca[dipole j]]
    # for a given scattering direction

    for j in eachindex(positions)
        jj = (j-1)*3+1:j*3    # find current dipole
        rj = positions[j]
        nrj = dot(n, rj)
        phase = exp(-1im * k * nrj)

        Esca += phase * G * dipoles[jj, :]
        # TODO Bsca as well

    end

    # Esca is now the net FF in given direction 
    return Esca
end




