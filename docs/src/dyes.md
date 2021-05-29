# Dimer of uniaxial dye molecules

This example considers the coupling between two uniaxial dye molecules in close proximity, with a side-by-side configuration.

```@example 1
using CoupledDipole
using StaticArrays
using DataFrames
using VegaLite

## materials
wavelength = collect(400:1:700.0)
media = Dict([("Rhodamine", alpha_bare), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

## dimer geometry
cl0 = cluster_single(0, 0, 1, 0, 0, 0, "Rhodamine", "point")
cl1 = cluster_dimer(0.8, 0, 0, 1, 0, 0,0, "Rhodamine", "point")

## incidence: along z, along x, along y
Incidence = [SVector(0,0,0), SVector(0,π/2,0), SVector(π/2,π/2,0)]

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

```

Note that there is essentially zero scattering since the molecules are _much_ smaller than the wavelength.