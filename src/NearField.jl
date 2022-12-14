"""
scattered_field(λ, probes, positions, dipoles)

- `λ`: wavelength
- `probes`: SVector of spatial positions where the field is to be evaluated
- `positions`: SVector of dipole positions
- `P`: 3Ndipx2Ninc matrix of self-consistent dipole moments

"""
function local_field(k, Rdip, Rpro, P, EinProbes)

    N_dip = length(Rdip)
    N_pro = length(Rpro)
    N_inc = size(P, 2)
    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / Z₀

    # Esca = zeros(eltype(P), (3N_pro, N_inc))
    # actually easier to work with array of vectors
    Esca = [@SVector(zeros(eltype(P), 3)) for _ ∈ 1:N_inc*N_pro]
    Bsca = similar(Esca)
    Etot = similar(Esca)
    Btot = similar(Esca)

    # for loop over N_probe locations
    for i = 1:N_pro
        ii = 3i-2:3i
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
                    k^2 * (I - nxn) -
                    (1 - im * k * rᵢⱼ) / (rᵢⱼ^2) * (I - 3 * nxn)
                )

            # EM 
            Bᵢⱼ = expikror * (nx - I) * (k^2 + im * k / rᵢⱼ)

            # contribution from j-th source dipole, for all incidences
            Etmp += Aᵢⱼ * P[jj, :]
            Btmp += Bᵢⱼ * P[jj, :]

        end
        # store scattered field 3-vector at probe location for each incidence
        # adding Ein is a bit awkward as it's stored in a 3NxNinc matrix like P
        for l ∈ 1:N_inc
            Esca[N_inc*(i-1)+l] = Etmp[:, l]
            Bsca[N_inc*(i-1)+l] = Z₀ * Btmp[:, l]
            Etot[N_inc*(i-1)+l] = Etmp[:, l] + EinProbes[ii, l]
            Btot[N_inc*(i-1)+l] = Btmp[:, l] + Y₀ * EinProbes[ii, l]
        end

    end

    return Esca, Bsca, Etot, Btot
end

ellipsoid(origin, size) = (origin[1] / size[1])^2 + (origin[2] / size[2])^2 + (origin[3] / size[3])^2

function is_inside(probe, positions, sizes, rotations)
    tests = map((p, s, r) -> ellipsoid(r * probe - p, s) <= 1, positions, sizes, rotations)
    return reduce(|, tests)
end


# internal_field(p, V, χ) = 1 / χ * p / V

function scattered_field(probe, k, Rdip, sizes, rotations, P)

    N_dip = length(Rdip)
    N_inc = size(P, 2)
    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / Z₀

    Esca = zeros(eltype(P), (3, N_inc))
    Bsca = zeros(eltype(P), (3, N_inc))

    inside = is_inside(probe, Rdip, sizes, rotations)
    # if inside
    #     Esca = internal_field(P[jj, :], prod(sizes[j]), χ)
    # else
    if !inside # we're outside particles
        for j = 1:N_dip
            jj = 3j-2:3j

            # source to probe
            rⱼ_to_rᵢ = probe - Rdip[j]
            rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
            n = rⱼ_to_rᵢ / rᵢⱼ
            nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n
            nx = SMatrix{3,3}(0, n[3], n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p

            expikror = exp(im * k * rᵢⱼ) / rᵢⱼ

            # EE 
            Aᵢⱼ =
                expikror * (
                    k^2 * (I - nxn) -
                    (1 - im * k * rᵢⱼ) / (rᵢⱼ^2) * (I - 3 * nxn)
                )

            # EM 
            Bᵢⱼ = expikror * (nx - I) * (k^2 + im * k / rᵢⱼ)

            # contribution from j-th source dipole, for all incidences
            Esca += Aᵢⱼ * P[jj, :]
            Bsca += Bᵢⱼ * P[jj, :]

        end
    end
    return Esca, Bsca, inside
end


function incident_field(Ejones, k, probe, IncidenceRotations)
    T = typeof(Ejones[1][1] * k * probe[1] * IncidenceRotations[1][1, 1])

    N_inc = length(IncidenceRotations)
    Ein = Array{T}(undef, (3, 2N_inc))

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

    for jj in eachindex(IncidenceRotations)
        Rm = IncidenceRotations[jj]
        k_hat = k * Rm[:, 3] # Rm * [0;0;1] == 3rd column
        E1_r = Rm * Evec1
        E2_r = Rm * Evec2
        kR = dot(k_hat, probe)
        Ein[:, jj] = E1_r * exp(im * kR)
        Ein[:, jj+N_inc] = E2_r * exp(im * kR)
    end

    return Ein
end

# x = -200.0:2.0:200
# probes = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]
# probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
function map_nf(probes,
    cl::Cluster,
    mat::Material,
    Incidence;
    polarisation="linear",
    prescription="kuwata")

    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / Z₀
    N_pro = length(probes)
    N_dip = length(cl.positions)
    N_inc = length(Incidence)

    ## low level stuff
    proto_r = cl.positions[1][1] # position type
    proto_a = RotMatrix(cl.rotations[1])[1, 1] # angle type
    proto_α = 0.1 + 0.1im # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α
    F = Matrix{T2}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{T2}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    # incident field
    if polarisation == "linear"
        Ejones = [
            SVector(1.0 + 0im, 0.0), # Jones vector, first polar
            SVector(0.0, 1.0 + 0im), # Jones vector, second polar
        ]
    elseif polarisation == "circular"
        Ejones = 1.0 / √2.0 .* [
            SVector(1im, 1.0), # Jones vector, first polar ↺
            SVector(1.0, 1im), # Jones vector, second polar ↻
        ]
    end

    #   solve the system for the polarisation
    ParticleRotations = map(RotMatrix, cl.rotations)
    IncidenceRotations = map(RotMatrix, Incidence)
    λ = mat.wavelengths[1]

    if length(mat.wavelengths) > 1
        @warn "this function expects a single wavelength; using $λ"
    end

    n_medium = mat.media["medium"](λ)
    k = n_medium * 2π / λ
    if cl.type == "point"
        Alpha = map(
            (m, s) ->
                alpha_scale(alpha_embed(mat.media[m](λ), n_medium), s),
            cl.materials,
            cl.sizes,
        )

    elseif cl.type == "particle"
        Epsilon = map(m -> mat.media[m](λ), cl.materials) # evaluate materials at wavelength
        Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ; prescription=prescription)
    end
    AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
    interaction_matrix_labframe!(F, k, cl.positions, AlphaBlocks)
    incident_field!(Ein, Ejones, k, cl.positions, IncidenceRotations)
    E = F \ Ein
    polarisation!(P, E, AlphaBlocks)

    # now the near-field part 

    T3 = typeof(proto_r)
    Esca = [@SMatrix(zeros(T2, 3, 2N_inc)) for _ ∈ 1:N_pro]
    Einc = similar(Esca)
    Bsca = similar(Esca)
    Etot = similar(Esca)
    Btot = similar(Esca)
    E² = Matrix{T3}(undef, N_pro, 2N_inc)
    B² = similar(E²)
    𝒞 = similar(E²)
    inside = Vector{Bool}(undef, N_pro)
    for i in eachindex(probes)

        Einc[i] = incident_field(Ejones, k, probes[i], IncidenceRotations)
        Esca[i], Bsca[i], inside[i] = scattered_field(probes[i], k, cl.positions, cl.sizes, ParticleRotations, P)
        Etot[i] = Einc[i] + Esca[i]
        Btot[i] = Y₀ * Einc[i] + Bsca[i]
        # scalar summaries, but for each incidence
        E²[i, :] = sum(abs2.(Etot[i]), dims=1)
        B²[i, :] = sum(abs2.(Btot[i]), dims=1)
        𝒞[i, :] = imag.(sum(conj.(Etot[i]) .* Btot[i], dims=1))

    end

    # for convenience, return the positions as a dataframe
    positions = DataFrame(reduce(vcat, transpose.(probes)), [:x, :y, :z])
    positions.inside .= inside

    return E², B², 𝒞, positions

end

