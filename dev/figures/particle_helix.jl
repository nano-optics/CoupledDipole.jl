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
wavelengths = collect(400:2:1000.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## dimer geometry
cl1 = cluster_helix(8, 20, 20, 40, 100, 400, π/2, 0)
cl2 = cluster_helix(8, 20, 20, 40, 100, 400, π/2, 0, "right")
cl0 = cluster_single(20, 20, 40)


oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl2, mat)

d0 = oa_df(oa0, mat.wavelengths)
d3 = oa_df(oa1, mat.wavelengths)
d4 = oa_df(oa2, mat.wavelengths)

d5 = [insertcols!(d3, :cluster => "helix", :hand => "left");
      insertcols!(d4, :cluster => "helix", :hand => "right");
      insertcols!(d0, :cluster => "single", :hand => "_")]


d5 |> @vlplot(
 width= 120,
 height =  100,
     mark = {:line},
     row = "type", column="variable",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q",
      strokeDash = "cluster:n",  color = "hand:n"}
 )



