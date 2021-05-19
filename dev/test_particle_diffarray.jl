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
wavelength = collect(450:2:850.0)
medium = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, medium)

## dimer geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl1 = CoupledDipole.cluster_array(900, 550, 50, 50, 60, 0,π/2,0)
cl2 = cluster_single(50, 50, 60, 0,π/2,0)


## incidence: along z
Incidence = [SVector(0,0,0)]

disp1 = spectrum_dispersion(cl1, mat, Incidence)
disp2 = spectrum_dispersion(cl2, mat, Incidence)

function dispersion_df(x, wavelength)

    e = insertcols!(
        DataFrame(x.extinction, :auto),
        :wavelength => wavelength,
        :crosstype => "extinction"
    )
    a = insertcols!(
        DataFrame(x.absorption, :auto),
        :wavelength => wavelength,
        :crosstype => "absorption"
    )

    s = insertcols!(
        DataFrame(x.scattering, :auto),
        :wavelength => wavelength,
        :crosstype => "scattering"
    )

    stack([e;a;s], Not([:wavelength,:crosstype]))

end

d1 = dispersion_df(disp1, mat.wavelength)
d2 = dispersion_df(disp2, mat.wavelength)

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
