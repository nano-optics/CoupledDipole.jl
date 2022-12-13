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

    return Esca, Bsca
end




function map_nf(probes,
    cl::Cluster,
    mat::Material,
    Incidence;
    polarisation="linear",
    prescription="kuwata")

    ## low level stuff
    proto_r = cl.positions[1][1] # position type
    proto_a = cl.rotations[1][1] # angle type
    proto_α = 0.1 + 0.1im # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α
    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)
    N_inc = length(Incidence)
    F = Matrix{T2}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{T2}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    Ejones = [SVector(1.0 + 0im, 0.0), SVector(0.0, 1.0 + 0im)]

    ParticleRotations = map(RotMatrix, cl.rotations) # now (active) Rotation objects
    IncidenceRotations = map(RotMatrix, Incidence) # now given as quaternions
    λ = mat.wavelengths[1]
    n_medium = mat.media["medium"](λ)
    kn = n_medium * 2π / λ
    Epsilon = map(m -> mat.media[m](λ), cl.materials)
    Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ)
    AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
    interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)
    incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)
    E = F \ Ein
    polarisation!(P, E, AlphaBlocks)

    x = -200.0:2.0:200
    probes = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]
    probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
    N_pro = length(probes)

    # when it was a matrix
    # Esca = scattered_field(kn, cl.positions, probes, P)
    # Isca = sum(reshape(abs2.(Esca), (3, N_pro * 2N_inc)), dims=1)
    # plot(collect(x), log10.(Isca[1:length(x)]))


    Esca, Bsca = scattered_field(kn, cl.positions, probes, P)

    E² = [sum(abs2.(E)) for E in Esca]
    B² = [sum(abs2.(B)) for B in Bsca]
    𝒞 = map((E, B) -> imag(E ⋅ B), Esca, Bsca)

    return E², B², 𝒞

end

