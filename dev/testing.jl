# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
# using Gadfly


# a =  cubature_sphere(12, "gl")


N_dip = 5
N_inc = 36
F = Matrix{Complex{Float64}}(I, 3N_dip, 3N_dip)
Ein = Array{Complex{Float64}}(undef, (3*N_dip, 2N_inc))
E = similar(Ein)
E = F\Ein


N_dip = 5
Ni = 36
R = [@SVector randn(3) for ii=1:N_dip]
Angles = [@SVector randn(3) for ii=1:N_dip]
kn = 0.2
Alpha = [@SVector randn(Complex{Float64},3) for ii=1:N_dip]

AlphaBlocks = [@SMatrix randn(Complex{Float64},3,3) for ii=1:N_dip]

# function alpha_blocks!(AlphaBlocks::Vector{SMatrix{3,3}},
#     Alpha::Vector{SVector{3}},
#     Angles::Vector{SVector{3}})


# alpha_blocks!(AlphaBlocks,Alpha,Angles)



A = Matrix{Complex{Float64}}(I, 3*N_dip, 3*N_dip)
propagator_freespace_labframe!(A,kn,R,AlphaBlocks)


# Evec = SVector{3, Complex{Float64}}([0.1; 0.1; 0.0])

Ejones = [SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0)]
Incidence = [@SVector randn(Float64,3) for ii=1:Ni]
IncidenceRotations = map(CoupledDipole.euler_active, Incidence)
Ein = Array{Complex{Float64}}(undef, (3*N_dip, 2Ni))
incident_field!(Ein, Ejones, kn, R, IncidenceRotations)
E = similar(Ein)
P = similar(Ein)
polarisation!(P, E, AlphaBlocks)

Cabs = Array{Float64}(undef, (2Ni))

absorption!(Cabs, kn, P, E)
extinction!(Cabs, kn, P, Ein)

cl = cluster_dimer(0.5, 1.0, 1.0, 1.0)

q = cubature_sphere(20)
scattering!(Cabs, cl.positions, q.nodes, q.weights, kn, P)


wavelength = collect(450:2:850.0)
epsilon = epsilon_Ag.(wavelength)

alpha_trace = vcat(alpha_bare.(wavelength)...)

# molecule
wavelength = collect(450:2:850.0)
medium = Dict([("alpha", alpha_bare), ("medium", x -> 1.33)])

mat = Material(wavelength, medium)
cl = cluster_dimer(0.5, 1.0, 1.0, 1.0)

# particle
wavelength = collect(450:2:850.0)
medium = Dict([("Au", CoupledDipole.epsilon_Au), ("medium", x -> 1.33)])

mat = Material(wavelength, medium)

cl = cluster_dimer(100, 20, 20, 40, π/4)
# cl = cluster_single(20, 20, 40)
#
# Alpha = alpha_spheroids(500, -10+1im, 1.33^2, cl.sizes)
#
# AlphaBlocks = [@SMatrix zeros(Complex{Real},3,3) for ii=1:N_dip]
#
# Ni=2
Incidence = [SVector(0,0,0),SVector(0,π/2,0)]


typeof(cl)
typeof(mat)
#
Juno.@enter spectrum_dispersion(cl, mat, Incidence)


testdisp = spectrum_dispersion(cl, mat, Incidence)

a = DataFrame(testdisp.scattering, :auto)
b = DataFrame(testdisp.absorption, :auto)
c = DataFrame(testdisp.extinction, :auto)
a[!,:wavelength] = mat.wavelength
b[!,:wavelength] = mat.wavelength
c[!,:wavelength] = mat.wavelength
da = stack(a, Not(:wavelength))
db = stack(b, Not(:wavelength))
dc = stack(c, Not(:wavelength))

@vlplot(
width= 400,
height =  300) +
@vlplot(data=da, #sca
    mark = {:line, strokeDash=(1,5)},
    encoding = {x = "wavelength:q", y = "value:q", color = "variable:n"}
) +
@vlplot(data=db, # abs
    mark = :line,
    encoding = {x = "wavelength:q", y = "value:q", color = "variable:n"}
) +
@vlplot(data=dc, # ext
    mark = {:line, strokeDash=(5,1)},
    encoding = {x = "wavelength:q", y = "value:q", color = "variable:n"}
)



d |> @vlplot(
    mark = :line,
    encoding = {x = "wavelength:q", y = "value", color = "variable:n"},
           tooltip={field="value", typ="quantitative"},
  width= 400,
  height =  300,
)


alpha = quadrature_lgwt( 5, 0, 2*π )

beta  = quadrature_lgwt( 4, 0, 1 )
a=2.0;b=1.0
nodes = hcat([SVector{3}(a,acos(2*b - 1),0.0) for b in beta.nodes, a in alpha.nodes]...)

weights = SVector([a*b for a in alpha.weights, b in beta.weights]...)

quad_inc = cubature_sphere(3, "gl")

cl = cluster_dimer(80, 20, 20, 40, pi/4)
testoa = spectrum_oa(cl, mat, "gl", 100)

clref = cluster_single(20,20,40,0,0,0)
# clref = cluster_single(20,20,40,0,π/2,0)

ref = spectrum_oa(clref, mat, "gl", 300)


a = DataFrame(dimer = testoa.average.extinction, single =  ref.average.extinction)
# a = DataFrame(dimer = testoa.dichroism.extinction, single =  ref.dichroism.extinction)
a[!,:wavelength] = mat.wavelength
d = stack(a, Not(:wavelength))

# plot(d, x=:wavelength, y=:value, colour=:variable, Geom.line)

d |> @vlplot(
    mark = :line,
    encoding = {x = "wavelength:q", y = "value", color = "variable:n"},
           tooltip={field="value", typ="quantitative"},
  width= 400,
  height = 300,
  autosize = "fit",
)
