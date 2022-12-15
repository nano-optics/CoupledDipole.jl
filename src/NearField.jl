
"""
scattered_field(probe, k, n, positions, sizes, rotations, P)

- `probe`: SVector, point where the scattered EM field is to be evaluated
- `k`: wavenumber
- `n_medium`: refractive index of medium 
- `positions`: SVector of dipole positions
- `sizes`: SVector of dipole sizes
- `rotations`: SVector of particle rotation matrices
- `P`: 3Ndipx2Ninc matrix of self-consistent dipole moments
- `evaluate_inside`: logical, whether to compute the (unphysical) scattered field inside particles

"""
function scattered_field(probe, k, n_medium, positions, sizes, rotations, P; evaluate_inside=true)

    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / 376.730313668 # H = Y₀ E
    c₀ = 299792458 # m/s
    ε₀ = 8.8541878128e-12

    N_dip = length(positions)
    N_inc = size(P, 2)

    Esca = zeros(eltype(P), (3, N_inc))
    Bsca = zeros(eltype(P), (3, N_inc))

    inside = is_inside(probe, positions, sizes, rotations)
    # if inside
    #     Esca = internal_field(P[jj, :], prod(sizes[j]), χ)
    # else
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
            Bᵢⱼ = expikror * (nx) * (k^2 + im * k / rᵢⱼ)

            # contribution from j-th source dipole, for all incidences
            Esca += Aᵢⱼ * P[jj, :]
            Bsca += Bᵢⱼ * P[jj, :]

        end
    end
    # B = Z ε₀ε Gm P̄ 
    scale = Z₀ * ε₀ * n_medium
    return Esca, scale * Bsca, inside
end


"""
incident_field(Ejones, k, probe, IncidenceRotations)

- `probe`: SVector, point where the scattered EM field is to be evaluated
- `k`: wavenumber
- n_medium
- `Ejones`: 2-Svectors defining 2 orthogonal Jones polarisations
- `IncidenceRotations`: `N_inc`-vector of rotation 3-Smatrices

"""
function incident_field(Ejones, k, n_medium, probe, IncidenceRotations)

    c₀ = 299792458 # m/s

    T = typeof(Ejones[1][1] * k * probe[1] * IncidenceRotations[1][1, 1])

    N_inc = length(IncidenceRotations)
    Ein = Array{T}(undef, (3, 2N_inc))
    Bin = similar(Ein)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector
    # B=(k x E)
    # Bx = - Ey   
    # By = Ex
    # Bz = 0 
    Bvec1 = SVector(-Ejones[1][2], Ejones[1][1], 0)
    Bvec2 = SVector(-Ejones[2][2], Ejones[2][1], 0)

    for jj in eachindex(IncidenceRotations)
        Rm = IncidenceRotations[jj]
        k_hat = k * Rm[:, 3] # Rm * [0;0;1] == 3rd column
        E1_r = Rm * Evec1
        E2_r = Rm * Evec2
        B1_r = Rm * Bvec1
        B2_r = Rm * Bvec2
        kR = dot(k_hat, probe)
        expikr = exp(im * kR)
        Ein[:, jj] = E1_r * expikr
        Ein[:, jj+N_inc] = E2_r * expikr
        Bin[:, jj] = n_medium / c₀ * B1_r * expikr
        Bin[:, jj+N_inc] = n_medium / c₀ * B2_r * expikr
    end

    return Ein, Bin
end

"""
map_nf(probes,
    cl::Cluster,
    mat::Material,
    Incidence;
    polarisation="linear",
    prescription="kuwata")

- `probes`: array of SVectors, points where the scattered EM field is to be evaluated
- `cl`: cluster of particles
- `mat`: dielectric functions
- `Incidence`: N_inc vector of quaternions describing incidence directions
- `polarisation`: incident field consists of 2 orthogonal "linear" or "circular" polarisations
- `prescription`: polarisability prescription for particles
- `evaluate_inside`: logical, whether to compute the (unphysical) scattered field inside particles

# example
# x = -200.0:2.0:200
# probes = SVector.(Iterators.product(x, x, zero(eltype(x))))[:]
# probes = SVector.(Iterators.product(x, zero(eltype(x)), zero(eltype(x))))[:]
"""
function map_nf(probes,
    cl::Cluster,
    mat::Material,
    Incidence;
    polarisation="linear",
    prescription="kuwata",
    evaluate_inside=true)

    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / 376.730313668 # H = Y₀ E
    c₀ = 299792458 # m/s

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
    Binc = similar(Esca)
    Bsca = similar(Esca)
    Etot = similar(Esca)
    Btot = similar(Esca)
    E² = Matrix{T3}(undef, N_pro, 2N_inc)
    B² = similar(E²)
    𝒞 = similar(E²)
    inside = Vector{Bool}(undef, N_pro)
    for i in eachindex(probes)

        Einc[i], Binc[i] = incident_field(Ejones, k, n_medium, probes[i], IncidenceRotations)
        Esca[i], Bsca[i], inside[i] = scattered_field(probes[i], k, n_medium, cl.positions, cl.sizes, ParticleRotations, P;
            evaluate_inside=evaluate_inside)
        Etot[i] = Einc[i] + Esca[i]
        Btot[i] = Binc[i] + Bsca[i]

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


ellipsoid(origin, size) = (origin[1] / size[1])^2 + (origin[2] / size[2])^2 + (origin[3] / size[3])^2

function is_inside(probe, positions, sizes, ParticleRotations)
    tests = map((p, s, r) -> ellipsoid(r' * (probe - p), s) <= 1, positions, sizes, ParticleRotations)
    return reduce(|, tests)
end

internal_field(p, V, χ) = 1 / χ * p / V
