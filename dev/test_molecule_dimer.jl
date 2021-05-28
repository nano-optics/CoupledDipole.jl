# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite

## materials
wavelength = collect(400:1:700.0)
media = Dict([("Rhodamine", alpha_bare), ("medium", x -> 1.33)])
mat = Material(wavelength, media)


## dimer geometry
cl0 = cluster_single(0, 0, 1, 0, 0, 0, "Rhodamine", "point")
cl1 = cluster_dimer(0.8, 0, 0, 1, 0, 0,0, "Rhodamine", "point")
cl2 = cluster_dimer(0.8, 0, 0, 1, π/4, 0,0, "Rhodamine", "point")

## incidence: along z, along x, along y
Incidence = [SVector(0,0,0),SVector(0,π/2,0),SVector(π/2,π/2,0)]

disp1 = spectrum_dispersion(cl0, mat, Incidence)
disp2 = spectrum_dispersion(cl1, mat, Incidence)

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



oa1 = spectrum_oa(cl0, mat)
oa2 = spectrum_oa(cl2, mat)


d3 = oa_df(oa1, mat.wavelengths)
d4 = oa_df(oa2, mat.wavelengths)

d5 = [insertcols!(d4, :cluster => "dimer");
     insertcols!(d3, :cluster => "single")]

d5 |> @vlplot(
 width= 400,
 height =  300,
     mark = {:line},
     row = "type",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )
