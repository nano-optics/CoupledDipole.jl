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

## dimer geometry
# cluster_dimer(d, a, b, c, ϕ=0.0, α_1=0.0, α_2=0.0, material="Au", type="particle")
cl1 = cluster_dimer(100.0, 20.0, 20.0, 35.0, π/4)
cl0 = cluster_single(20.0, 20.0, 35.0)

## incidence: along z, along x, along y
angles = range(0,π/2,7)
Incidence = QuatRotation.([RotX.(angles);RotY.(angles);RotZ.(angles)])

disp1 = spectrum_dispersion(cl0, mat, Incidence, polarisations="circular")
disp2 = spectrum_dispersion(cl1, mat, Incidence, polarisations="circular")

dd1 = dispersion_df(disp1, mat.wavelengths, format="wide")
dd2 = dispersion_df(disp2, mat.wavelengths, format="wide")

oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)

d0 = oa_df(oa0, mat.wavelengths)
d3 = oa_df(oa1, mat.wavelengths)

# d5 = [insertcols!(dd1, :cluster => "dimer", :hand => "left");
#       insertcols!(dd2, :cluster => "single", :hand => "_")]


using DataFramesMeta
dich = @chain dd2 begin
    @rsubset :crosstype == "extinction"
    @transform :dichroism = :pol1 .- :pol2
end

rotations = DataFrame(:angle => repeat(angles*180/π,outer=3), 
                      :axis => repeat(["x","y","z"], inner=length(angles)),
                      :angle_id => eachindex(Incidence))

toplot = leftjoin(dich, rotations, on=:angle_id)

toplot |> @vlplot(
 width= 120,
 height =  100,
     mark = {:line},
     row = "axis", 
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "dichroism:q",
      strokeDash = "cluster:n",  color = "angle:n"}
 )
