# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

# using LinearAlgebra
using StaticArrays
# using FastGaussQuadrature
# using DataFrames
using VegaLite

## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

## dimer geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl1 = CoupledDipole.cluster_array(900, 550, 50, 50, 60, 0,π/2,0)
cl2 = cluster_single(50, 50, 60, 0,π/2,0)


## incidence: along z
Incidence = [SVector(0,0,0)]

disp1 = spectrum_dispersion(cl1, mat, Incidence)
disp2 = spectrum_dispersion(cl2, mat, Incidence)


d1 = dispersion_df(disp1, mat.wavelengths)
d2 = dispersion_df(disp2, mat.wavelengths)

d = [insertcols!(d1, :cluster => "dimer");
     insertcols!(d2, :cluster => "single")]

 @vlplot(data=d,
 width= 400,
 height =  300,
     mark = {:line},
     row = "crosstype",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )
