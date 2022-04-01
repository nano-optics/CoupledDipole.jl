# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using CoupledDipole
using StaticArrays
using DataFrames
using VegaLite

## materials
wavelengths = collect(400:2:1000.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## array geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl0 = cluster_single(80, 40, 40)
cl1 = cluster_array(200, 550, 80, 40, 40)

## incidence: along z
Incidence = [SVector(0,0,0)]

disp1 = spectrum_dispersion(cl0, mat, Incidence)
disp2 = spectrum_dispersion(cl1, mat, Incidence)


d1 = dispersion_df(disp1, mat.wavelengths)
d2 = dispersion_df(disp2, mat.wavelengths)

d = [insertcols!(d1, :cluster => "single");
     insertcols!(d2, :cluster => "array")]

 @vlplot(data=d,
 width= 400,
 height =  300,
     mark = {:line},
     row = "crosstype",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )
