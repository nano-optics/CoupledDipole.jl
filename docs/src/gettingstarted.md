# Getting Started


## Installation

```julia
Pkg.add("CoupledDipole")
```

## Example

A typical simulation requires two inputs: a Cluster, describing the configuration of particles, and a Material, describing the wavelength-dependent dielectric functions of all media.

```@example 1
using CoupledDipole
using DataFrames
using VegaLite

## cluster geometry
cl1 = cluster_helix(4, 20, 20, 40, 200, 400, π/4, 0)
cl0 = cluster_single(20, 20, 40) # reference: single-particle

## materials
wavelengths = collect(400:2:1000.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

```

From these two objects we can simply call a high-level function to simulate optical properties. The following lines simulate the orientation-averaged optical response,


```@example 1
oa0 = spectrum_oa(cl0, mat) # reference: just one particle
oa1 = spectrum_oa(cl1, mat)
```

From there we combine the cross-sections in long-format dataframes for plotting,

```@example 1
d0 = oa_df(oa0, mat.wavelengths)
d1 = oa_df(oa1, mat.wavelengths)

d2 = [insertcols!(d1, :cluster => "helix");
      insertcols!(d0, :cluster => "single")]

d2 |> @vlplot(
 width= 120,
 height =  100,
     mark = {:line, clip = true},
     row = "type", column="variable",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q",
      strokeDash = "cluster:n",  color = "hand:n"}
 )

```


![](helix-plot.svg)
