
## high-level wrappers

"""
    incident_field_pw!(Ein,
        Ejones,
        kn, R::Vector{SVector{3}},
        IncidenceRotations)

Incident field at particle positions

- `Ein`: `3N_dip x N_inc` matrix, right-hand side of coupled-dipole system
- `Ejones`: tuple of 2 2-Svectors defining 2 orthogonal Jones polarisations
- `kn`: wavenumber in incident medium
- `R`: `N_dip`-vector of 3-Svectors of particle positions
- `IncidenceRotations`: `N_inc`-vector of rotation 3-Smatrices

"""
function incident_field_pw!(Ein, Ejones, kn, R, IncidenceRotations)

    N_inc = length(IncidenceRotations)
    N_dip = length(R)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector

    for jj in eachindex(IncidenceRotations)
        Rm = IncidenceRotations[jj]
        # Rm * [0;0;1] == 3rd column
        k_hat = kn * Rm[:, 3]
        E1_r = Rm * Evec1
        E2_r = Rm * Evec2
        for kk = 1:N_dip

            expikr = exp(im * dot(k_hat, R[kk]))

            # use views to avoid copies (?)
            # Ev1 = @view Ein[kk*3-2:kk*3, jj]
            # Ev2 = @view Ein[kk*3-2:kk*3, jj+N_inc]
            # fill!(Ev1, E1_r * exp(im * kR))
            # fill!(Ev2, E2_r * exp(im * kR))

            @views Ein[kk*3-2:kk*3, jj] = E1_r * expikr
            @views Ein[kk*3-2:kk*3, jj+N_inc] = E2_r * expikr
        end
    end
    return Ein
end

function incident_field_pw(Ejones, kn, R, IncidenceRotations)

    N_dip = length(R)
    N_inc = length(IncidenceRotations)

    # what is the type of arrays to initialise?
    proto_r = R[1][1] # position type
    proto_E = Ejones[1][1] # E type
    proto_a = RotMatrix(IncidenceRotations[1])[1, 1] # angle type
    T = typeof(im * kn * proto_a * proto_E * proto_r) # Einc ~ Ejones . exp(ikr) 

    Ein = Array{T}(undef, (3N_dip, 2N_inc))
    incident_field_pw!(Ein, Ejones, kn, R, IncidenceRotations)
    return Ein
end


"""
incident_field_probe(Ejones, k, probe, IncidenceRotations)

evaluates the E and B fields at a probe location, for mapping purposes
Note that much of the functionality is common with incident_field_pw

- `probe`: SVector, point where the scattered EM field is to be evaluated
- `k`: wavenumber
- n_medium
- `Ejones`: 2-Svectors defining 2 orthogonal Jones polarisations
- `IncidenceRotations`: `N_inc`-vector of rotation 3-Smatrices

"""
function incident_field_probe(Ejones, k, n_medium, probe, IncidenceRotations)

    # TODO refactor to avoid duplication with incident_field_pw

    c₀ = 299792458 # m/s

    T = typeof(Ejones[1][1] * k * probe[1] * IncidenceRotations[1][1, 1])

    N_inc = length(IncidenceRotations)
    Ein = Array{T}(undef, (3, 2N_inc))
    Bin = similar(Ein)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector
    # B⃗ = k⃗ × E⃗
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
    incident_field_pw!(Ein, source::SVector{3}, kn, positions::Vector{SVector{3}})

Incident field at particle positions, from a dipole source

- `Ein`: `3N_dip x N_inc` matrix, right-hand side of coupled-dipole system
- `source`: 3-Svector, source position
- `kn`: wavenumber in incident medium
- `positions`: `N_dip`-vector of 3-Svectors of particle positions

"""
function incident_field_dip!(Ein, source, k, positions)

    N_dip = length(positions)

    # consider 3 unit dipoles for each location
    p₁ = SVector(1.0, 0.0, 0.0) # x
    p₂ = SVector(0.0, 1.0, 0.0) # y
    p₃ = SVector(0.0, 0.0, 1.0) # z


    for j in eachindex(positions)
        jj = 3j-2:3j # corresponding rows in Ein

        # source to particle
        rⱼ_to_rᵢ = positions[j] - source
        rᵢⱼ = norm(rⱼ_to_rᵢ, 2)
        n = rⱼ_to_rᵢ / rᵢⱼ
        nxn = n * transpose(n) # (n ⊗ n) p = (n⋅p) n

        kr = k * rᵢⱼ
        k2expikror = k^2 * exp(im * kr) / rᵢⱼ

        Aᵢⱼ = k2expikror * (
            (I - nxn) + (im / kr - 1 / kr^2) * (I - 3 * nxn)
        )
        # EM magnetic field unused for now
        # nx = SMatrix{3,3}(0, n[3], -n[2], -n[3], 0, n[1], n[2], -n[1], 0) # n × p
        # Bᵢⱼ = expikror * nx * (k^2 + im * k / rᵢⱼ)

        # source incident on j-th particle, for 3 directions of source

        @views Ein[jj, 1] = Aᵢⱼ * p₁
        @views Ein[jj, 2] = Aᵢⱼ * p₂
        @views Ein[jj, 3] = Aᵢⱼ * p₃
    end
    return Ein
end

