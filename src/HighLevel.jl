
## high-level wrappers

"""
     spectrum_dispersion(cl::Cluster, mat::Material,
                         Incidence, N_sca::Int=36)

Simulating far-field cross-sections for multiple wavelengths and directions of incidence

- `cl`: cluster of particles
- `mat`: dielectric functions
- `Incidence`: N_inc vector of 3-Svectors of incidence Euler angles
- `N_sca`: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_dispersion(
    cl::Cluster,
    mat::Material,
    Incidence;
    N_sca::Int = 36,
    method = "direct",
)

    # what is the type of arrays to initialise?
    proto_r = cl.positions[1][1]
    proto_k = 2π / mat.wavelengths[1]
    T = typeof(proto_k * proto_r)
    # note: cross-sections are typeof(imag(P*E)), which boils down to T, hopefully

    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)
    N_inc = length(Incidence)

    # scattering angles for Csca
    quad_sca = cubature_sphere(N_sca, "gl")

    # initialise
    # https://discourse.julialang.org/t/initialise-array-without-specific-types/60204/3
    F = Matrix{Complex{T}}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{Complex{T}}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    # AlphaBlocks = [@SMatrix zeros(Complex{T}, 3, 3) for ii = 1:N_dip]

    # incident field
    Ejones = [
        SVector{2}(1.0 + 0im, 0.0), # Jones vector, first polar
        SVector{2}(0.0, 1.0 + 0im), # Jones vector, second polar
    ]

    # store all rotation matrices
    ParticleRotations = map(euler_passive, cl.angles)
    IncidenceRotations = map(euler_active, Incidence)
    ScatteringVectors = map(euler_unitvector, quad_sca.nodes)
    # NOTE: we only need the third column to rotate kz, should specialise

    ## loop over wavelengths

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Array{T}(undef, (N_lam, 2 * N_inc))
    cabs = similar(cext)
    csca = similar(cext)
    tmpcext = Vector{T}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)
    tmpcsca = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelengths[ii]
        n_medium = mat.media["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            α_name = cl.material # e.g. "alpha" to refer to mat Dict
            α_bare = mat.media[α_name](λ)
            α = alpha_embedded(α_bare, n_medium)
            Alpha = alpha_rescale_molecule(α, cl.sizes)

        elseif cl.type == "particle"

            ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
            ε = mat.media[ε_name](λ)
            Alpha = alpha_spheroids(λ, ε, n_medium^2, cl.sizes)

        end

        # update the rotated blocks
        # Remark: if we cared about memory more than speed
        # we'd compute rotation matrices on the fly
        # alpha_blocks!(AlphaBlocks, Alpha, cl.angles)
        # but instead we've prestored them, since wavelength-independent
        AlphaBlocks =
            map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)

        # Interaction matrix (F = I - G0 alpha_eff)
        propagator_freespace_labframe!(F, kn, cl.positions, AlphaBlocks)

        # update the incident field
        # Remark: similarly, pre-computed rotations
        # incident_field!(Ein, Ejones, kn, cl.positions, Incidence)
        # but instead we've prestored them, since wavelength-independent
        incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)


        # solve

        if method == "direct" # brute-force solve linear system

            E = F \ Ein
            polarisation!(P, E, AlphaBlocks)

            # cross-section for cubature angles
            extinction!(tmpcext, kn, P, Ein)

        elseif method == "oos" # udpate E, P, α_ext iteratively by order-of-scattering
            # 
            # E[:] = Ein
            # polarisation!(P, E, AlphaBlocks)
            # extinction!(tmpcext, kn, P, Ein)
            iterate_field!(
                E,
                P,
                tmpcext,
                Ein,
                F,
                kn,
                AlphaBlocks
            )
        end

        # remaining cross-sections for cubature angles
        absorption!(tmpcabs, kn, P, E)
        scattering!(
            tmpcsca,
            cl.positions,
            ScatteringVectors,
            quad_sca.weights,
            kn,
            P,
        )

        cext[ii, :] = tmpcext
        cabs[ii, :] = tmpcabs
        csca[ii, :] = tmpcsca
    end

    CrossSections(1 / N_dip * cext, 1 / N_dip * cabs, 1 / N_dip * csca)

end


"""
     spectrum_oa(cl::Cluster, mat::Material,
                Cubature = "gl", N_inc::Int = 36, N_sca::Int=36)

Orientation-averaged far-field cross-sections for multiple wavelengths

- `cl`: cluster of particles
- `mat`: dielectric functions
- `Cubature`: spherical cubature method
- `N_inc`: number of incident angles for spherical cubature
- `N_sca`: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_oa(
    cl::Cluster,
    mat::Material;
    cubature = "gl",
    N_inc = 36,
    N_sca = 36,
    method = "direct"
)

    quad_inc = cubature_sphere(N_inc, cubature)
    quad_sca = cubature_sphere(N_sca, cubature)

    # setting up constants

    # what is the type of arrays to initialise?
    proto_r = cl.positions[1][1]
    proto_k = 2π / mat.wavelengths[1]
    T = typeof(proto_k * proto_r)
    # note: cross-sections are typeof(imag(P*E)), which boils down to T, hopefully

    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)
    N_inc = length(quad_inc.weights) # update with actual number of angles

    F = Matrix{Complex{T}}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{Complex{T}}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    AlphaBlocks = [@SMatrix zeros(Complex{T}, 3, 3) for ii = 1:N_dip]
    # block_array = BlockArray{Float32}(undef_blocks, [1,2], [3,2])
    # setblock!(block_array, rand(2,2), 2, 1)
    # block_array[Block(1, 1)]

    # incident field
    Ejones = 1.0 / sqrt(2.0) .* [
        SVector{2}(1im, 1.0), # Jones vector, first polar
        SVector{2}(1.0, 1im), # Jones vector, second polar
    ]

    # store all rotation matrices
    ParticleRotations = map(euler_passive, cl.angles)
    IncidenceRotations = map(euler_active, quad_inc.nodes)
    ScatteringVectors = map(euler_unitvector, quad_sca.nodes)

    # average both polarisations, so divide by two
    weights1 = 0.5 * vcat(quad_inc.weights, quad_inc.weights) #  standard cross sections
    weights2 = vcat(quad_inc.weights, -quad_inc.weights) #  dichroism

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Vector{T}(undef, N_lam)
    cabs = similar(cext)
    csca = similar(cext)
    dext = similar(cext)
    dabs = similar(cext)
    dsca = similar(cext)

    tmpcext = Vector{T}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)
    tmpcsca = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelengths[ii]
        n_medium = mat.media["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            α_name = cl.material # e.g. "alpha" to refer to mat Dict
            α_bare = mat.media[α_name](λ)
            α = alpha_embedded(α_bare, n_medium)
            Alpha = alpha_rescale_molecule(α, cl.sizes)

        elseif cl.type == "particle"

            ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
            ε = mat.media[ε_name](λ)
            Alpha = alpha_spheroids(λ, ε, n_medium^2, cl.sizes)

        end

        # update the rotated blocks
        # alpha_blocks!(AlphaBlocks, Alpha, cl.angles)
        AlphaBlocks =
            map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)

        propagator_freespace_labframe!(F, kn, cl.positions, AlphaBlocks)

        # incident_field!(Ein, Ejones, kn, cl.positions, quad_inc.nodes)
        incident_field!(Ein, Ejones, kn, cl.positions, IncidenceRotations)

        # solve

        if method == "direct" # brute-force solve linear system

            E = F \ Ein
            polarisation!(P, E, AlphaBlocks)

            # cross-section for cubature angles
            extinction!(tmpcext, kn, P, Ein)

        elseif method == "oos" # udpate E, P, α_ext iteratively by order-of-scattering

            E = Ein
            polarisation!(P, E, AlphaBlocks)
            extinction!(tmpcext, kn, P, Ein)
            iterate_field!(
                E,
                P,
                tmpcext,
                Ein,
                F,
                kn,
                AlphaBlocks
            )

        end

        # remaining cross-sections for cubature angles
        absorption!(tmpcabs, kn, P, E)
        scattering!(
            tmpcsca,
            cl.positions,
            ScatteringVectors,
            quad_sca.weights,
            kn,
            P,
        )


        #  perform cubature for angular averaging
        cext[ii] = dot(tmpcext, weights1)
        cabs[ii] = dot(tmpcabs, weights1)
        csca[ii] = dot(tmpcsca, weights1)
        dext[ii] = dot(tmpcext, weights2)
        dabs[ii] = dot(tmpcabs, weights2)
        dsca[ii] = dot(tmpcsca, weights2)
    end

    # csca = cext - cabs
    # dsca = dext - dabs

    (
        average = CrossSections(
            1 / N_dip * cext,
            1 / N_dip * cabs,
            1 / N_dip * csca,
        ),
        dichroism = CrossSections(
            1 / N_dip * dext,
            1 / N_dip * dabs,
            1 / N_dip * dsca,
        ),
    )


end
