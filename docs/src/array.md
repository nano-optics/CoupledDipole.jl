# Diffractive array of gold nanorods

This example considers a square array of 200 Au nanorods with a large pitch (550nm); this configuration can lead to a diffractive coupling effect that strongly modifies the resonance of individual rods, despite their wide separation.

```@example 1
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

```

We'll do the simulation at normal incidence, and simulate the properties of a single rod for reference,

```@example 1

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

```

Note that scattering cross-sections are terribly inaccurate here; because of the large size of the cluster the numerical quadrature would require a large number of scattering angles. It would be preferable to evaluate scattering as the difference between extinction and absorption.


