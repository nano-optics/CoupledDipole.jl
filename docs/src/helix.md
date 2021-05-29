# Helix of gold nanorods

This example considers a short helix of gold nanorods, and simulates its orientation-averaged circular dichroism. Two helices are modelled, with opposite handedness, to verify that they produce mirror-image dichroism spectra.


```@example 1
using CoupledDipole
using DataFrames
using VegaLite

## materials
wavelengths = collect(400:2:1000.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## cluster geometry
cl1 = cluster_helix(8, 20, 20, 40, 100, 400, π/2, 0)
cl2 = cluster_helix(8, 20, 20, 40, 100, 400, π/2, 0, "right")
cl0 = cluster_single(20, 20, 40)


oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl2, mat)

d0 = oa_df(oa0, mat.wavelengths)
d1 = oa_df(oa1, mat.wavelengths)
d2 = oa_df(oa2, mat.wavelengths)

d = [insertcols!(d1, :cluster => "helix", :hand => "left");
      insertcols!(d2, :cluster => "helix", :hand => "right");
      insertcols!(d0, :cluster => "single", :hand => "_")]


d |> @vlplot(
 width= 120,
 height =  100,
     mark = {:line},
     row = "type", column="variable",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q",
      strokeDash = "cluster:n",  color = "hand:n"}
 )

```

