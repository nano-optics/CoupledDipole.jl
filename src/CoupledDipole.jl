module CoupledDipole


using LinearAlgebra
using StaticArrays
using BlockArrays
using Base.Threads
using FastGaussQuadrature
using Gadfly
using DataFrames

# learnings
# Types: avoid specifying, so that the code works with more generic types e.g. dual numbers needed by AD
# Svectors, Sarrays: best for efficiency. Initialise them with a tuple, not an array, to avoid needless conversion
# Block matrix: to keep its initialisation type-free, probably want to use BlockArray package, defining blocks in a tuple
# function defs: don't give types, instead suggest them in docstrings for reference of intended usage
# broadcasting: easier to construct most functions with scalars, and broadcast on array arguments when needed
#

struct Cluster{T1,T2,T3}
    positions::Vector{SVector{3,T1}}
    angles::Vector{SVector{3,T2}}
    sizes::Vector{SVector{3,T3}}
end

struct Material{T1,T2,T3}
    wavelength::Vector{T1}
    alpha::Vector{T2}
    medium::T3
end

struct CrossSections{T}
    extinction::Array{T}
    absorption::Array{T}
    scattering::Array{T}
end

# struct Cubature{T}
#     nodes::SMatrix{T}
#     weights::SVector{T}
# end

## cluster definitions



@doc raw"""
    cluster_single(a::T, b::T, c::T, α::T = 0.0, β::T = 0.0, γ::T = 0.0) where T <: Real

Particle cluster consisting of a single particle at the origin
- a,b,c: semi-axis along x,y,z
- α,β,γ: Euler angles

@example
   cluster_single(1.0,2.0,3.0)
"""
function cluster_single(a, b, c, α = 0.0, β = 0.0, γ = 0.0)
    sizes = [SVector{3}(a, b, c)]
    positions = [SVector{3}(0.0, 0.0, 0.0)]
    angles = [SVector{3}(α, β, γ)]
    Cluster(positions, angles, sizes)
end


"""
    cluster_dimer(d::T, a::T, b::T, c::T, dihedral::T = 0.0, α_1::T = 0.0, α_2::T = 0.0) where T <: Real

Particle cluster consisting of 2 identical particles separated along z
- a,b,c: semi-axis along x,y,z
- dihedral: angle between both particles seen along the z-axis
- α_1,α_2: polar Euler angle of each particle

@example
   cluster_dimer(10,1,2,3)
"""
function cluster_dimer(d, a, b, c, dihedral = 0.0, α_1 = 0.0, α_2 = 0.0)
    sizes = [SVector{3}(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector{3}(0.0, 0.0, z) for z in (d/2, -d/2)]
    angles = [SVector{3}(dihedral, α_1, 0.0),
              SVector{3}(0.0, α_2, 0.0)]
    Cluster(positions, angles, sizes)
end


## material functions

"""
    epsilon_Ag(λ::T) where T <: Real

Drude model for the dielectric function of silver in the visible region
- λ: wavelength in nm

@example
   epsilon_Ag(632.8)
"""
function epsilon_Ag(λ)
    4*(1.0 - 1.0 / (282.0^2 * (1 / λ^2 + im / (17000 * λ))))
end

"""
    lorentzian(λ::T, α_k::T, λ_k::T, µ_k::T) where T <: Real

Complex Lorentz function, to describe complex polarisabilities
- λ: wavelength in nm
- α_k: oscillator strength in S.I. units
- λ_k: oscillator wavelength in nm
- µ_k: damping in S.I. units

@example
   lorentzian(632.8)
"""
function lorentzian(λ, α_k = 5.76e-38, λ_k = 526.0, µ_k = 1.0e4)
    - α_k*λ_k/µ_k * (1.0 - 1.0 / ( 1.0 - (λ_k/ λ)^2 - 1im*(λ_k^2 /(µ_k*λ))))
end

"""
    alpha_bare(λ::T, α_∞::T, α_k::Array{T}, λ_k::Array{T}, µ_k::Array{T}) where T <: Real

Complex scalar polarisability, as sum of lorentz oscillators
- λ: wavelength in nm
- α_k: oscillator strength(s) in S.I. units
- λ_k: oscillator wavelength(s) in nm
- µ_k: damping(s) in S.I. units

Default values mimic the main resonance of Rhodamine 700

@example
   alpha_bare(632.8)
"""
function alpha_bare(λ, α_∞=9.6e-39,
    α_k=5.76e-38, λ_k=526.0, µ_k=1.0e4)

    ε_0 = 8.8541878128e-12
    nm3 = 1e27

    α = α_∞
    for kk in 1:length(α_k)
        α += lorentzian(λ, α_k[kk], λ_k[kk], µ_k[kk])
    end
    prefact = nm3/(4.0*pi*ε_0)
    prefact * α
end


"""
    alpha_embedded(α::Complex{T}, medium::T) where T <: Real

Complex scalar polarisability, as sum of lorentz oscillators
- α: bare polarisabilty
- medium: refractive index of embedding medium

Default values mimic the main resonance of Rhodamine 700

@example
   alpha_embedded(alpha_bare(632.8))
"""
function alpha_embedded(α, medium = 1.33)
    ε_m = medium^2
    L =  (ε_m + 2) / 3
    1/ε_m * L^2 * α
end


"""
    alpha_rescale(alpha, sizes::Vector{SVector{3}})

Principal polarisability components of a particle, rescaled along each principal axis
- α: scalar polarisabilty
- sizes: array of 3-vectors to scale along each principal axis

@example
   sizes = [SVector{3}(1.0, 2.0, 3.0) for i in 1:4]
   alpha_rescale(alpha_bare(632.8), sizes)
"""
function alpha_rescale(alpha, sizes)
    alpha .* (sizes ./ sum.(sizes))
end



## orientation-averaging


@doc raw"""
    quadrature_lgwt(N::Int, a::Real, b::Real)

N-point Gauss-Legendre quadrature over [a,b] interval
- N: number of nodes
- a,b: bounds

$$
\int _{a}^{b}f(x)\,dx={\frac {b-a}{2}}\int _{-1}^{1}f\left({\frac {b-a}{2}}x+{\frac {a+b}{2}}\right)\,dx.
$$

@example
   quadrature_lgwt(6,0,3)
"""
function quadrature_lgwt(N, a, b)
    n,w = FastGaussQuadrature.gausslegendre(N)
    nn = (b-a)/2.0 * n .+ (a+b) / 2.0
    ww = (b-a)/2.0 * w
    (nodes = nn, weights = ww)
end

"""
    cubature_sphere(N::Int, method::String)

N-point Gauss-Legendre quadrature over [a,b] interval
- N: number of nodes
- method: cubature method (only 'gl' currently implemented)

Returns a Cubature object containing 2 arrays (Nx3 nodes and Nx1 weights)
Note: using array instead of tuple for weights because we'll probably use them in a scalar product in orientation-averaging
For nodes there is less of a reason, but it can be convenient to visualise the nodes.

The cubature is normalised by 4π such that a unit integrand approximates 1.

@example
   cubature_sphere(36*18)
"""
function cubature_sphere(N, method="gl")
    #might have slightly more than N total points
    rndN = Integer(ceil(sqrt(N/2.0)))

    alpha = quadrature_lgwt( 2*rndN, 0, 2*pi )
    beta  = quadrature_lgwt( rndN, 0, 1 )

    nodes = [SVector{3}(a,acos(2*b - 1),0.0) for b in beta.nodes, a in alpha.nodes]
    weights = hcat([a*b for a in alpha.weights, b in beta.weights]...)

    (nodes = nodes, weights = 1/(2*pi)*weights)
end

## utils

"""
    euler_active(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- φ: Euler angle (longitude, in [0,2*pi])
- θ: Euler angle (colatitude, in [0,2*pi])
- ψ: Euler angle (rotation around z", in [0,2*pi])


@example
   euler_active(pi/2,0,0)
"""
function euler_active(φ, θ, ψ)

    cosφ = cos(φ); cosψ = cos(ψ); cosθ = cos(θ)
    sinφ = sin(φ); sinψ = sin(ψ); sinθ = sin(θ)

    R = SMatrix{3,3}(
    cosφ*cosθ*cosψ - sinφ*sinψ, -cosφ*cosθ*sinψ - sinφ*cosψ, cosφ*sinθ,
    sinφ*cosθ*cosψ + cosφ*sinψ, -sinφ*cosθ*sinψ + cosφ*cosψ, sinφ*sinθ,
    -sinθ*cosψ, sinθ*sinψ, cosθ)

    return R

end

"""
    euler_passive(φ::Real, θ::Real, ψ::Real)

3D rotation matrix

- φ: Euler angle (longitude, in [0,2*pi])
- θ: Euler angle (colatitude, in [0,2*pi])
- ψ: Euler angle (rotation around z", in [0,2*pi])


@example
   euler_passive(pi/2,0,0)
"""
function euler_passive(φ, θ, ψ)

    cosφ = cos(φ); cosψ = cos(ψ); cosθ = cos(θ)
    sinφ = sin(φ); sinψ = sin(ψ); sinθ = sin(θ)

    R = SMatrix{3,3}(
     cosφ*cosθ*cosψ - sinφ*sinψ,  sinφ*cosθ*cosψ + cosφ*sinψ, -sinθ*cosψ,
    -cosφ*cosθ*sinψ - sinφ*cosψ, -sinφ*cosθ*sinψ + cosφ*cosψ,  sinθ*sinψ,
     cosφ*sinθ, sinφ*sinθ, cosθ)

    return R

end



## core functions



"""
    alpha_blocks(AlphaBlocks::Vector{SMatrix{3,3}},
                 Alpha::Vector{SVector{3}},
                 Angles::Vector{SVector{3}})

Polarisability blocks for the interaction matrix

- AlphaBlocks: N-vector of 3x3 Smatrices (polarisability tensors in the lab frame)
- Alpha: N-vector of principal polarisability components for each particle
- Angles: Euler angles for rotating each particle in the lab frame

"""
function alpha_blocks!(AlphaBlocks, Alpha, Angles)

    N = length(Angles)

    # loop over N particles
    for ii in 1:N
        Rm = euler_passive(Angles[ii]...)
        AlphaBlocks[ii] =  Rm' * (Alpha[ii] .* Rm)
    end

    return AlphaBlocks
end

"""
    propagator_freespace_labframe(A,
        kn, R::Vector{SVector{3}},
        AlphaBlocks::Vector{SMatrix{3,3}})

Polarisability blocks for the interaction matrix

- A: 3Nx3N interaction matrix
- kn: wavenumber in incident medium
- R: N-vector of 3-Svectors of particle positions
- AlphaBlocks: N-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

# try BlockArray.jl
# Construct a BlockArray from blocks
# mortar((A,B),(C,D))

"""
function propagator_freespace_labframe!(A, kn, R, AlphaBlocks)

    N = length(R)

    # nested for loop over N dipoles
    for jj in 1:N

        ind_jj = 3jj-2:3jj

        for kk in (jj+1):N

            ind_kk = 3kk-2:3kk

            rk_to_rj = R[jj] - R[kk]
            rjk = norm(rk_to_rj, 2)
            rjkhat = rk_to_rj / rjk
            rjkrjk =  rjkhat * transpose(rjkhat)

            Ajk = exp(im*kn*rjk) / rjk * (kn*kn*(rjkrjk - I) +
            (im*kn*rjk - 1.0) / (rjk*rjk) * (3*rjkrjk - I))

            # assign blocks
            A[ind_jj, ind_kk] = Ajk * AlphaBlocks[kk]
            A[ind_kk, ind_jj] = transpose(Ajk) * AlphaBlocks[jj]

        end
    end

    return A
end



"""
    incident_field(Ein,
        Ejones,
        kn, R::Vector{SVector{3}},
        Incidence)

Incident field at particle positions

- Ein: 3NxNi matrix, right-hand side of coupled-dipole system
- Ejones: tuple of 2 2-Svectors defining 2 orthogonal Jones polarisations
- kn: wavenumber in incident medium
- R: N-vector of 3-Svectors of particle positions
- Incidence: Ni-vector of 3-Svectors of incidence angles

# try BlockArray.jl
# Construct a BlockArray from blocks
# mortar((A,B),(C,D))

"""
function incident_field!(Ein, Ejones, kn, R, Incidence)

    Ni = length(Incidence)
    N = length(R)

    Evec1 = SVector(Ejones[1][1], Ejones[1][2], 0) # 3-vector
    Evec2 = SVector(Ejones[2][1], Ejones[2][2], 0) # 3-vector
    Evec_r = similar(Evec1)

    for jj in eachindex(Incidence)

        Rm = euler_active(Incidence[jj]...)
        # Rot * [0.0;0.0;1.0] == 3d column
        k_hat = kn * transpose(Rm[:,3])
        for kk in 1:N
            Ein[kk*3-2:kk*3,jj]    = (Rm * (Evec1 * exp(im * k_hat * R[kk])))
            Ein[kk*3-2:kk*3,jj+Ni] = (Rm * (Evec2 * exp(im * k_hat * R[kk])))
        end
    end
    return Ein
end


"""
    polarisation!(P, E, AlphaBlocks)

Incident field at particle positions

- P: 3NxNi matrix, polarisations for all incidences
- E: 3NxNi matrix, total field for all incidences
- AlphaBlocks: N-vector of 3x3 Smatrices (polarisability tensors in the lab frame)

"""
function polarisation!(P, E, AlphaBlocks)

    # loop over N particles
    for ii in eachindex(AlphaBlocks)
        ind = 3ii-2:3ii
        P[ind,:] = AlphaBlocks[ii] * E[ind,:] # all incidence angles at once
    end

    return P
end


## cross-sections

"""
    extinction(kn::Real, P::Array{Complex}, Ein::Array{Complex})

Extinction cross-section for each incident angle

- kn: wavenumber in incident medium
- P:   3NxNi matrix, polarisations for all incidences
- Ein: 3NxNi matrix, incident field for all incidences

"""
function extinction!(cext, kn, P, Ein)

    Ni = size(P,2)
    # cext = Vector{Float64}(undef, Ni)

    for jj in 1:Ni
        cext[jj] = 4*pi*kn*imag(dot(Ein[:,jj], P[:,jj])) # E^* P
    end

    return  cext
end


"""
    absorption(kn::Real, P::Array{Complex}, E::Array{Complex})

Absorption cross-section for each incident angle

- kn: wavenumber in incident medium
- P:   3NxNi matrix, polarisations for all incidences
- E: 3NxNi matrix, total field for all incidences

"""
function absorption!(cabs, kn, P, E)

    Ni = size(P, 2)
    # cabs = Vector{Float64}(undef, Ni)

    for jj in 1:Ni
        cabs[jj] = 4*pi*kn*imag(dot(E[:,jj], P[:,jj]) -
        kn^3 * 2/3 * real(dot(P[:,jj],P[:,jj])))
    end

    return  cabs
end


# note: need also Csca as numerical cubature, for consistency check with Ext-Abs

## high-level wrappers

"""
     spectrum_dispersion(cl::Cluster, material::Material,
                         Incidence, N_sca::Int=36)

Simulating far-field cross-sections for multiple wavelengths and directions of incidence

- cl: cluster of particles
- material: dielectric functions
- Incidence: Ni vector of 3-Svectors of incidence Euler angles
- N_sca: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_dispersion(cl::Cluster, material::Material,
    Incidence, N_sca::Int=36)

    T = typeof(cl.positions[1][1]) # type used for array inits below
    N_dip = length(cl.positions)
    N_lambda = length(material.wavelength)
    N_inc = length(Incidence)

# initialise

    # https://discourse.julialang.org/t/initialise-array-without-specific-types/60204/3
    F = Matrix{Complex{T}}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{Complex{T}}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    AlphaBlocks = [@SMatrix zeros(Complex{T},3,3) for ii=1:N_dip]

    # incident field
    kn = material.medium .* 2*pi./material.wavelength
    Ejones =  [SVector{2}(1.0 + 0im, 0.0), # Jones vector, first polar
               SVector{2}(0.0, 1.0 + 0im)] # Jones vector, second polar

    # scattering angles for Csca2
    # quad_sca = cubature_sphere(N_sca, 'gl')

    ## loop over wavelengths

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Array{T}(undef,(2*N_inc, N_lambda))
    cabs = similar(cext)

    for ii in 1:N_lambda

        # TODO switch between dyes (point dipoles) and particles
        # TODO allow dispersive medium

        alpha = alpha_embedded(material.alpha[ii], material.medium)
        Alpha = alpha_rescale(alpha, cl.sizes)
        alpha_blocks!(AlphaBlocks, Alpha, cl.angles)

        # Interaction matrix (A = I - G0 alpha_eff)
        propagator_freespace_labframe!(F, kn[ii], cl.positions, AlphaBlocks)

        incident_field!(Ein, Ejones, kn[ii], cl.positions, Incidence)

        # @show size(Ein)
        # @show size(F)
        # solve
        E = F \ Ein
        polarisation!(P, E, AlphaBlocks)

        # cross-sections for multiple angles
        extinction!(cext[:,ii], kn[ii], P, Ein)
        absorption!(cabs[:,ii], kn[ii], P, E)
    end

    CrossSections(1/N_dip*cext, 1/N_dip*cabs, 1/N_dip*(cext-cabs))

end


"""
     spectrum_oa(cl::Cluster, material::Material,
                Cubature = "gl", N_inc::Int = 36, N_sca::Int=36)

Orientation-averaged far-field cross-sections for multiple wavelengths

- cl: cluster of particles
- material: dielectric functions
- Cubature: spherical cubature method
- N_inc: number of incident angles for spherical cubature
- N_sca: number of scattering angles for spherical cubature estimate of σ_sca

"""
function spectrum_oa(cl::Cluster, material::Material,
                     cubature = "gl", N_inc = 36, N_sca = 36)


    quad_inc = cubature_sphere(N_inc, cubature)
    quad_sca = cubature_sphere(N_sca, cubature)

    # setting up constants
    T = typeof(cl.positions[1][1]) # type used for array inits below
    N_dip = length(cl.positions)
    N_lambda = length(material.wavelength)
    N_inc = length(quad_inc.weights) # update with actual number of angles

    F = Matrix{Complex{T}}(I, 3N_dip, 3N_dip) # type inferred from cl.positions
    Ein = Array{Complex{T}}(undef, (3N_dip, 2N_inc))
    E = similar(Ein)
    P = similar(Ein)

    AlphaBlocks = [@SMatrix zeros(Complex{T},3,3) for ii=1:N_dip]
    # block_array = BlockArray{Float32}(undef_blocks, [1,2], [3,2])
    # setblock!(block_array, rand(2,2), 2, 1)
    # block_array[Block(1, 1)]

    # incident field
    kn =  material.medium .* 2*pi ./ material.wavelength

    # incident field
    kn = material.medium .* 2*pi./material.wavelength
    Ejones =  1.0/sqrt(2.0).*[SVector{2}(1im,1.0), # Jones vector, first polar
                              SVector{2}(1.0,1im)] # Jones vector, second polar

    # average both polarisations, so divide by two
    weights1 = 0.5*vcat(quad_inc.weights,  quad_inc.weights) #  standard cross sections
    weights2 =     vcat(quad_inc.weights, -quad_inc.weights) #  dichroism

    # use type T for containers, inferred from positions
    # this is to be compatible with autodiff when optimising cluster geometry
    cext = Vector{T}(undef,N_lambda)
    cabs = similar(cext)
    csca2 = similar(cext)
    dext = similar(cext)
    dabs = similar(cext)
    dsca2 = similar(cext)

    tmpcext = Vector{T}(undef,2N_inc)
    tmpcabs = similar(tmpcext)

    for ii in 1:N_lambda

        # TODO switch between dyes and particles
        # TODO allow dispersive medium
        alpha = alpha_embedded(material.alpha[ii], material.medium)
        Alpha = alpha_rescale(alpha, cl.sizes)
        alpha_blocks!(AlphaBlocks, Alpha, cl.angles)

        propagator_freespace_labframe!(F, kn[ii], cl.positions, AlphaBlocks)

        incident_field!(Ein, Ejones, kn[ii], cl.positions, quad_inc.nodes)

        E = F \ Ein
        polarisation!(P, E, AlphaBlocks)

        # cross-sections for Cubature angles
        extinction!(tmpcext,kn[ii], P, Ein)
        absorption!(tmpcext,kn[ii], P, E)

        #  perform Cubature for angular averaging
        cext[ii]  = dot(tmpcext,  weights1)
        cabs[ii]  = dot(tmpcabs,  weights1)
        # csca2[ii] = dot(tmpsca2, weights1)
        dext[ii]  = dot(tmpcext,  weights2)
        dabs[ii]  = dot(tmpcabs,  weights2)
        # dsca2[ii] = dot(tmpsca2, weights2)
    end

    csca = cext - cabs
    dsca = dext - dabs

    #  return a matrix
    #  also normalise by number of dipoles
    CrossSections(1/N_dip*cext, 1/N_dip*cabs, 1/N_dip*(cext-cabs))

end


end
