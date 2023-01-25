
## high-level wrappers

"""
     spectrum_dispersion(cl::Cluster, mat::Material,
                         Incidence, N_sca::Int=36)

Simulating far-field cross-sections for multiple wavelengths and directions of incidence

- `cl`: cluster of particles
- `mat`: dielectric functions
- `Incidence`: N_inc vector of quaternions describing incidence directions
- `polarisation`: incident field consists of 2 orthogonal "linear" or "circular" polarisations
- `N_sca`: number of scattering angles for spherical cubature estimate of σ_sca
- `prescription`: polarisability prescription for particles
- `method`: direct or iterative solver

"""
function spectrum_dispersion(
    cl::Cluster,
    mat::Material,
    Incidence;
    N_sca::Int=36,
    polarisation="linear",
    prescription="kuwata",
    method="direct"
)

    # what is the type of arrays to initialise?
    proto_r = cl.positions[1][1] # position type
    # CHECK: would picking the first element of quaternion be a sensible idea
    proto_a = RotMatrix(cl.rotations[1])[1, 1] # angle type
    proto_α = proto_r * (1.0 + 1.0im) # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    # note: cross-sections are typeof(imag(P*E)), which boils down to T1, hopefully
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α

    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)
    N_inc = length(Incidence)

    # scattering angles for Csca
    quad_sca = cubature_sphere(N_sca, "gl")

    # initialise
    # https://discourse.julialang.org/t/initialise-array-without-specific-types/60204/3
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

    # store all rotation matrices
    # ParticleRotations = map(euler_passive, cl.angles)
    ParticleRotations = map(RotMatrix, cl.rotations) # now (active) Rotation objects
    # IncidenceRotations = map(euler_active, Incidence) # old Euler, no longer needed
    IncidenceRotations = map(RotMatrix, Incidence) # now given as quaternions
    ScatteringVectors = map(euler_unitvector, quad_sca.nodes)
    # TODO we only need the third column to rotate kz, should specialise

    ## loop over wavelengths

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Array{T1}(undef, (N_lam, 2N_inc))
    cabs = similar(cext)
    csca = similar(cext)
    tmpcext = Vector{T1}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)
    tmpcsca = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelengths[ii]
        n_medium = mat.media["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            Alpha = map(
                (m, s) ->
                    alpha_scale(alpha_embed(mat.media[m](λ), n_medium), s),
                cl.materials,
                cl.sizes,
            )

        elseif cl.type == "particle"

            # Alpha = map(
            #     (m, s) -> alpha_kuwata(λ, mat.media[m](λ), n_medium^2, s),
            #     cl.materials,
            #     cl.sizes,
            # )

            Epsilon = map(m -> mat.media[m](λ), cl.materials) # evaluate materials at wavelength
            Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ; prescription=prescription)
        end

        #  R is passive
        AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
        # when R is active 
        # AlphaBlocks = map((R, A) -> R * (diagm(A) * R'), ParticleRotations, Alpha)

        # Interaction matrix (F = I - G0 alpha_eff)
        interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)

        # update the incident field
        # Remark: similarly, pre-computed rotations
        # incident_field_pw!(Ein, Ejones, kn, cl.positions, Incidence)
        # but instead we've prestored them, since wavelength-independent
        incident_field_pw!(Ein, Ejones, kn, cl.positions, IncidenceRotations)


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
            iterate_field!(E, P, tmpcext, Ein, F, kn, AlphaBlocks)
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
- `prescription`: polarisability prescription for particles
- `method`: solution of linear system (direct or iterative)

"""
function spectrum_oa(
    cl::Cluster,
    mat::Material;
    cubature="gl",
    N_inc=36,
    N_sca=36,
    prescription="kuwata",
    method="direct"
)

    quad_inc = cubature_sphere(N_inc, cubature)
    quad_sca = cubature_sphere(N_sca, cubature)

    # setting up constants

    # what is the type of arrays to initialise?
    proto_r = cl.positions[1][1] # position type
    proto_a = cl.rotations[1][1] # angle type
    proto_α = 0.1 + 0.1im # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    # note: cross-sections are typeof(imag(P*E)), which boils down to T1, hopefully
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α


    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)
    N_inc = length(quad_inc.weights) # update with actual number of angles

    F = Matrix{T2}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{T2}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    # AlphaBlocks = [@SMatrix zeros(T2, 3, 3) for ii = 1:N_dip]
    # block_array = BlockArray{Float32}(undef_blocks, [1,2], [3,2])
    # setblock!(block_array, rand(2,2), 2, 1)
    # block_array[Block(1, 1)]

    # incident field
    Ejones = 1.0 / √2.0 .* [
        SVector(1im, 1.0), # Jones vector, first polar ↺
        SVector(1.0, 1im), # Jones vector, second polar ↻
    ]

    # store all rotation matrices
    # ParticleRotations = map(euler_passive, cl.angles) # old version
    ParticleRotations = map(RotMatrix, cl.rotations) # now (active) Rotation objects
    # IncidenceRotations = map(euler_active, quad_inc.nodes) # old version, replaced by RotZYZ
    IncidenceRotations = map(RotZYZ, quad_inc.nodes)
    ScatteringVectors = map(euler_unitvector, quad_sca.nodes)

    # average both polarisations, so divide by two
    weights1 = 0.5 * vcat(quad_inc.weights, quad_inc.weights) #  standard cross sections
    weights2 = vcat(quad_inc.weights, -quad_inc.weights) #  dichroism

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Vector{T1}(undef, N_lam)
    cabs = similar(cext)
    csca = similar(cext)
    dext = similar(cext)
    dabs = similar(cext)
    dsca = similar(cext)

    tmpcext = Vector{T1}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)
    tmpcsca = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelengths[ii]
        n_medium = mat.media["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            # old version: only one material per cluster
            # α_name = cl.material # e.g. "alpha" to refer to mat Dict
            # α_bare = mat.media[α_name](λ)
            # α = alpha_embedded(α_bare, n_medium)
            # Alpha = alpha_rescale_molecule(α, cl.sizes)

            Alpha = map(
                (m, s) ->
                    alpha_scale(alpha_embed(mat.media[m](λ), n_medium), s),
                cl.materials,
                cl.sizes,
            )

        elseif cl.type == "particle"

            # ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
            # ε = mat.media[ε_name](λ)
            # Alpha = alpha_particles(λ, ε, n_medium^2, cl.sizes)
            # Alpha = map(
            #     (m, s) -> alpha_kuwata(λ, mat.media[m](λ), n_medium^2, s),
            #     cl.materials,
            #     cl.sizes,
            # )
            Epsilon = map(m -> mat.media[m](λ), cl.materials) # evaluate materials at wavelength
            Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ; prescription=prescription)
        end

        #  R is passive
        AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)
        # when R is active 
        # AlphaBlocks = map((R, A) -> R * (diagm(A) * R'), ParticleRotations, Alpha)

        interaction_matrix_labframe!(F, kn, cl.positions, AlphaBlocks)

        # incident_field_pw!(Ein, Ejones, kn, cl.positions, quad_inc.nodes)
        incident_field_pw!(Ein, Ejones, kn, cl.positions, IncidenceRotations)

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
            iterate_field!(E, P, tmpcext, Ein, F, kn, AlphaBlocks)

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


        # cabs2[] = oa_c_abs(k, G0, invA, diagInvAlpha)

    end


    (
        average=CrossSections(
            1 / N_dip * cext,
            1 / N_dip * cabs,
            1 / N_dip * csca,
        ),
        dichroism=CrossSections(
            1 / N_dip * dext,
            1 / N_dip * dabs,
            1 / N_dip * dsca,
        ),
    )


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
- `prescription`: polarisability prescription for particles
- `method`: solution of linear system (direct or iterative)

"""
function spectrum_oa_analytical(
    cl::Cluster,
    mat::Material;
    prescription="kuwata",
    method="direct"
)

    # setting up constants

    # what is the type of arrays to initialise?
    proto_r = cl.positions[1][1] # position type
    proto_a = cl.rotations[1][1] # angle type
    proto_α = 0.1 + 0.1im # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    # note: cross-sections are typeof(imag(P*E)), which boils down to T1, hopefully
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α

    N_dip = length(cl.positions)
    N_lam = length(mat.wavelengths)

    ParticleRotations = map(RotMatrix, cl.rotations)

    cext = Vector{T1}(undef, N_lam)
    cabs = similar(cext)
    csca = similar(cext)

    for i = 1:N_lam
        λ = mat.wavelengths[i]
        n_medium = mat.media["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            Alpha = map(
                (m, s) ->
                    alpha_scale(alpha_embed(mat.media[m](λ), n_medium), s),
                cl.materials,
                cl.sizes,
            )

        elseif cl.type == "particle"

            Epsilon = map(m -> mat.media[m](λ), cl.materials)
            Alpha = alpha_particles(Epsilon, cl.sizes, n_medium^2, λ; prescription=prescription)

        end

        AlphaBlocks = map((R, A) -> R' * (diagm(A) * R), ParticleRotations, Alpha)

        tmpcabs, tmpcext = oa_analytical(kn, cl.positions, AlphaBlocks)

        cabs[i] = tmpcabs
        cext[i] = tmpcext

    end

    csca = cext - cabs

    # unimplemented
    # miss = Vector{Union{Float64,Missing}}(missing, N_lam)
    (
        average=CrossSections(
            1 / N_dip * cext,
            1 / N_dip * cabs,
            1 / N_dip * csca,
        ),
        dichroism=CrossSections( # TODO
            0 * cext,
            0 * cext,
            0 * cext,
        ),
    )


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
    evaluate_inside=true,
    return_fields=false)

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
    incident_field_pw!(Ein, Ejones, k, cl.positions, IncidenceRotations)
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
        # TODO refactor to use same incident field function as PW, varargout E,B
        Einc[i], Binc[i] = incident_field_probe(Ejones, k, n_medium, probes[i], IncidenceRotations)
        Esca[i], Bsca[i], inside[i], id = scattered_field(probes[i], k, n_medium, cl.positions, cl.sizes, ParticleRotations, P, Epsilon;
            evaluate_inside=evaluate_inside)
        if inside[i] # internal field deduced from P is the total field
            Etot[i] = Esca[i]
        else
            Etot[i] = Einc[i] + Esca[i]
        end
        Btot[i] = Binc[i] + Bsca[i]

        # scalar summaries, but for each incidence
        E²[i, :] = sum(abs2.(Etot[i]), dims=1)
        B²[i, :] = sum(abs2.(Btot[i]), dims=1)
        𝒞[i, :] = c₀ / n_medium * imag.(sum(conj.(Etot[i]) .* Btot[i], dims=1))

    end

    # for convenience, return the positions as a dataframe
    positions = DataFrame(reduce(vcat, transpose.(probes)), [:x, :y, :z])
    positions.inside .= inside
    if return_fields
        return Einc, Binc, Esca, Bsca, Etot, Btot, positions
    end
    return E², B², 𝒞, positions

end



"""
scattering_pattern(probes,
    cl::Cluster,
    mat::Material,
    source;
    polarisation="linear",
    prescription="kuwata")

- `directions`: array of SVectors, directions where the scattered EM field is to be evaluated
- `cl`: cluster of particles
- `mat`: dielectric functions
- `source`: SVector of dipole source position
- `prescription`: polarisability prescription for particles

# example
# φ = range(0,2π,36)
# θ = range(0,π,18)
# directions = SVector.(Iterators.product(φ, θ, 0.0))[:]
"""
function scattering_pattern(directions,
    cl::Cluster,
    mat::Material,
    source;
    prescription="kuwata")

    Z₀ = 376.730313668 # free-space impedance
    Y₀ = 1 / 376.730313668 # H = Y₀ E
    c₀ = 299792458 # m/s

    N_dir = length(directions)
    N_dip = length(cl.positions)


    ## low level stuff
    proto_r = cl.positions[1][1] # position type
    proto_a = RotMatrix(cl.rotations[1])[1, 1] # angle type
    proto_α = 0.1 + 0.1im # dummy complex polarisability
    proto_k = 2π / mat.wavelengths[1]
    T1 = typeof(proto_k * proto_r * imag(proto_α * proto_a)) #
    T2 = typeof(proto_k * proto_r + proto_α * proto_a) # blocks are ~ exp(ikr) or R * α
    F = Matrix{T2}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{T2}(undef, (3N_dip, 3)) # 3 orthogonal dipoles considered
    E = similar(Ein)
    P = similar(Ein)


    #   solve the system for the polarisation
    ParticleRotations = map(RotMatrix, cl.rotations)
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

    # incident field from dipole sources
    incident_field_dip!(Ein, source, k, cl.positions)

    E = F \ Ein
    polarisation!(P, E, AlphaBlocks)

    # now the scattered field pattern

    T3 = typeof(proto_r)
    E² = Matrix{T3}(undef, N_dir, 3)
    p₁ = SVector(1.0, 0.0, 0.0) # x
    p₂ = SVector(0.0, 1.0, 0.0) # y
    p₃ = SVector(0.0, 0.0, 1.0) # z
    Psources = hcat(p₁, p₂, p₃)
    dipoles = vcat(Psources, P)
    # @info display(dipoles)
    positions = copy(cl.positions)
    insert!(positions, 1, source)
    for i in eachindex(directions)

        Esca = far_field(directions[i], k, positions, dipoles)

        # scalar summaries, but for each emission direction
        E²[i, :] = sum(abs2.(Esca), dims=1)

    end

    return E²

end
