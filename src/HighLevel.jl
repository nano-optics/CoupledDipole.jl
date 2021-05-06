
## high-level wrappers

"""
     spectrum_dispersion(cl::Cluster, material::Material,
                         Incidence, N_sca::Int=36)

Simulating far-field cross-sections for multiple wavelengths and directions of incidence

- cl: cluster of particles
- material: dielectric functions
- Incidence: N_inc vector of 3-Svectors of incidence Euler angles
- N_sca: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_dispersion(
    cl::Cluster,
    mat::Material,
    Incidence,
    N_sca::Int = 36,
)

    T = typeof(cl.positions[1][1]) # type used for array inits below
    N_dip = length(cl.positions)
    N_lam = length(mat.wavelength)
    N_inc = length(Incidence)

    # scattering angles for Csca
    quad_sca = cubature_sphere(N_sca, cubature)

    # initialise
    # https://discourse.julialang.org/t/initialise-array-without-specific-types/60204/3
    F = Matrix{Complex{T}}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{Complex{T}}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    AlphaBlocks = [@SMatrix zeros(Complex{T}, 3, 3) for ii = 1:N_dip]

    # incident field
    Ejones = [
        SVector{2}(1.0 + 0im, 0.0), # Jones vector, first polar
        SVector{2}(0.0, 1.0 + 0im), # Jones vector, second polar
    ]


    ## loop over wavelengths

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Array{T}(undef, (N_lam, 2 * N_inc))
    cabs = similar(cext)
    tmpcext = Vector{T}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelength[ii]
        n_medium = mat.medium["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            α_name = cl.material # e.g. "alpha" to refer to mat Dict
            α_bare = mat.medium[α_name](λ)
            α = alpha_embedded(α_bare, n_medium)
            Alpha = alpha_rescale_molecule(α, cl.sizes)

        elseif cl.type == "particle"

            ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
            ε = mat.medium[ε_name](λ)
            Alpha = alpha_spheroids(λ, ε, n_medium^2, cl.sizes)

        end

        # update the rotated blocks
        # TODO the rotation matrices should be computed only once and stored
        # since wavelength-independent
        # (if we don't care about memory)
        alpha_blocks!(AlphaBlocks, Alpha, cl.angles)

        # Interaction matrix (A = I - G0 alpha_eff)
        propagator_freespace_labframe!(F, kn, cl.positions, AlphaBlocks)

        incident_field!(Ein, Ejones, kn, cl.positions, Incidence)

        # solve
        # TODO provide iterative solver alternative
        E = F \ Ein
        polarisation!(P, E, AlphaBlocks)

        # cross-sections for multiple angles
        extinction!(tmpcext, kn, P, Ein)
        absorption!(tmpcabs, kn, P, E)
        cext[ii, :] = tmpcext
        cabs[ii, :] = tmpcabs
    end

    CrossSections(1 / N_dip * cext, 1 / N_dip * cabs, 1 / N_dip * (cext - cabs))

end


"""
     spectrum_oa(cl::Cluster, mat::Material,
                Cubature = "gl", N_inc::Int = 36, N_sca::Int=36)

Orientation-averaged far-field cross-sections for multiple wavelengths

- cl: cluster of particles
- mat: dielectric functions
- Cubature: spherical cubature method
- N_inc: number of incident angles for spherical cubature
- N_sca: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_oa(
    cl::Cluster,
    mat::Material,
    cubature = "gl",
    N_inc = 300,
    N_sca = 36,
)

    quad_inc = cubature_sphere(N_inc, cubature)
    quad_sca = cubature_sphere(N_sca, cubature)

    # setting up constants
    T = typeof(cl.positions[1][1]) # type used for array inits below
    N_dip = length(cl.positions)
    N_lam = length(mat.wavelength)
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

    # average both polarisations, so divide by two
    weights1 = 0.5 * vcat(quad_inc.weights, quad_inc.weights) #  standard cross sections
    weights2 = vcat(quad_inc.weights, -quad_inc.weights) #  dichroism

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Vector{T}(undef, N_lam)
    cabs = similar(cext)
    csca2 = similar(cext)
    dext = similar(cext)
    dabs = similar(cext)
    dsca2 = similar(cext)

    tmpcext = Vector{T}(undef, 2N_inc)
    tmpcabs = similar(tmpcext)

    for ii = 1:N_lam
        λ = mat.wavelength[ii]
        n_medium = mat.medium["medium"](λ)
        kn = n_medium * 2π / λ

        if cl.type == "point"

            α_name = cl.material # e.g. "alpha" to refer to mat Dict
            α_bare = mat.medium[α_name](λ)
            α = alpha_embedded(α_bare, n_medium)
            Alpha = alpha_rescale_molecule(α, cl.sizes)

        elseif cl.type == "particle"

            ε_name = cl.material # e.g. "Au" to refer to epsilon_Au in mat Dict
            ε = mat.medium[ε_name](λ)
            Alpha = alpha_spheroids(λ, ε, n_medium^2, cl.sizes)

        end

        # update the rotated blocks
        # TODO the rotation matrices should be computed only once and stored
        # since wavelength-independent
        # (if we don't care about memory)
        alpha_blocks!(AlphaBlocks, Alpha, cl.angles)

        propagator_freespace_labframe!(F, kn, cl.positions, AlphaBlocks)

        incident_field!(Ein, Ejones, kn, cl.positions, quad_inc.nodes)

        # solve
        # TODO provide iterative solver alternative
        E = F \ Ein
        polarisation!(P, E, AlphaBlocks)

        # cross-sections for cubature angles
        extinction!(tmpcext, kn, P, Ein)
        absorption!(tmpcext, kn, P, E)

        #  perform cubature for angular averaging
        cext[ii] = dot(tmpcext, weights1)
        cabs[ii] = dot(tmpcabs, weights1)
        # csca2[ii] = dot(tmpsca2, weights1)
        dext[ii] = dot(tmpcext, weights2)
        dabs[ii] = dot(tmpcabs, weights2)
        # dsca2[ii] = dot(tmpsca2, weights2)
    end

    csca = cext - cabs
    dsca = dext - dabs

    (
        average = CrossSections(
            1 / N_dip * cext,
            1 / N_dip * cabs,
            1 / N_dip * (cext - cabs),
        ),
        dichroism = CrossSections(
            1 / N_dip * dext,
            1 / N_dip * dabs,
            1 / N_dip * (dext - dabs),
        ),
    )


end
