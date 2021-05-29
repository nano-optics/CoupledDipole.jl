# Chiral dimer of gold nanorods

This example considers the "fingers-crossed" configuration of two gold nanorods arranged in a chiral dimer, with a `pi/4` dihedral angle. The dipole-dipole interaction leads to circular dichroism with a characteristic bisignate lineshape.
 

```@example 1
using CoupledDipole
using DataFrames
using VegaLite


## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

## dimer geometry
cl0 = cluster_single(20, 20, 40)
cl1 = cluster_dimer(80, 20, 20, 40, π/4)

```

The following lines simulate the orientation-averaged optical response,


```@example 1
oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)

```

From there we combine the cross-sections in long-format dataframes for plotting,

```@example 1

d0 = oa_df(oa0, mat.wavelengths)
d1 = oa_df(oa1, mat.wavelengths)

d = [insertcols!(d1, :cluster => "dimer");
     insertcols!(d0, :cluster => "single")]

d5 |> @vlplot(
 width= 400,
 height =  300,
     mark = {:line},
     row = "type",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )

```

