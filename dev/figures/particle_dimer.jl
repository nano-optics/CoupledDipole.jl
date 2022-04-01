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
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

## reference geometry
cl0 = cluster_single(20, 20, 40)

function model(;d=80, ϕ=0) 
    cl = cluster_dimer(d, 20, 20, 40, ϕ*π/180)
    oa = spectrum_oa(cl, mat)
    oa_df(oa, mat.wavelengths)
end

# model()

params = expand_grid(d = range(80, 200, step = 20), ϕ = (0,45))

all = pmap_df(params, p -> model(;p...))

s = spectrum_oa(cl0, mat)
single = oa_df(s, mat.wavelengths)

# using AlgebraOfGraphics, CairoMakie
# set_aog_theme!()


# xy = data(all) * mapping(:wavelength, :value, color = :d => nonnumeric, col=:variable, row=:type, linestyle=:ϕ  => nonnumeric)
# layer = visual(Lines)
# AlgebraOfGraphics.draw(layer * xy, facet = (;  linkyaxes = :none))

x = repeat(LinRange(-5,5,100), inner=2)
df = DataFrame(x=x, y=x.^2, z=x.^3, col=repeat(1:2, outer=100))

df |>
@vlplot(facet={column={field=:col, typ=:nominal}}) +
  (
    @vlplot(x=:x) +
    @vlplot(:line, y=:y) +
    @vlplot(mark = {:line,color=:black}, strokeDash = {value = [2,2]}, y=:z)
)

single |> 
@vlplot(facet={column={field=:variable, typ=:nominal},
row={field=:type, typ=:nominal}}) +
(
    @vlplot() +
    @vlplot(mark = {:line},data=all,row=:type, column=:variable,
     encoding = {x = "wavelength:q", y = "value:q", color = "d:n", strokeDash="ϕ:n"}
     ) +
 @vlplot(mark = {:line,color=:black},data=single,
 strokeDash = {value = [2,2]},row=:type, column=:variable,
encoding = {x = "wavelength:q", y = "value:q"}
   )
)

@vlplot(mark = {:line},data=all, row=:type, column=:variable,
encoding = {x = "wavelength:q", y = "value:q", color = "d:n", strokeDash="ϕ:n"}
)


using VegaDatasets
dataset("weather.csv") |>
@vlplot(repeat={column=[:temp_max,:precipitation,:wind]}) +
(
    @vlplot() +
    @vlplot(
        :line,
        y={field={repeat=:column},aggregate=:mean,typ=:quantitative},
        x="month(date):o",
        detail="year(date):t",
        color=:location,
        opacity={value=0.2}
    ) +
    @vlplot(
        :line,
        y={field={repeat=:column},aggregate=:mean,typ=:quantitative},
        x="month(date):o",
        color=:location
    )
)

# using Gadfly

# plot(all, x=:wavelength, y=:value, linestyle = :ϕ, color=:d,
#     ygroup=:type, xgroup=:variable,
#      Geom.subplot_grid(Geom.line, free_y_axis=true))


using RDatasets: dataset
using AlgebraOfGraphics, CairoMakie
mpg = dataset("ggplot2", "mpg");
cols = mapping(:Displ, :Hwy);
grp = mapping(color = :Cyl => nonnumeric);
scat = visual(Scatter)
pipeline = cols * scat
data(mpg) * pipeline |> draw