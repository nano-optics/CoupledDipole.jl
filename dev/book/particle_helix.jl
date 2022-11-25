# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using Rotations

## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

# cl <- cluster_helix(5, R0=20, pitch=30,
#                     delta=pi/2, delta0=0, right=TRUE,
#                     a=15/2, b=15/2, c=15,
#                     angles="helix")

# nn <- 7
# Incidence <- rep(seq(0,pi/2,length=nn), 3)


## dimer geometry
# cluster_helix(N, a, b, c, R, Λ, δ, δ_0=0, handedness="left",
#     material="Au", type="particle")
# d, a, b, c, ϕ=0.0
cl1 <- cluster_dimer(100, 20.0, 20.0, 35.0, π/4)
cl1 = cluster_helix(5, 15.0/2, 15.0/2, 15.0, 20, 30, π/2, 0)
cl2 = cluster_helix(5, 15.0/2, 15.0/2, 15.0, 20, 30, π/2, 0, "right")
cl0 = cluster_single(15.0/2, 15.0/2, 15.0)

## incidence: along z, along x, along y
angles = range(0,π/2,7)
Incidence = QuatRotation.([RotX.(angles);RotY.(angles);RotZ.(angles)])

disp1 = spectrum_dispersion(cl0, mat, Incidence)
disp2 = spectrum_dispersion(cl1, mat, Incidence)

dd1 = dispersion_df(disp1, mat.wavelengths)
dd2 = dispersion_df(disp2, mat.wavelengths)

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
     row = "type", column="crosstype",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q",
      strokeDash = "cluster:n",  color = "hand:n"}
 )
