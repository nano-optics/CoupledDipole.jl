# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using Rotations
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


# dimer_model <- function(d, orientation = c("head-to-tail", "side-by-side"), 
#                         material, ...){

#   pb$tick()

#   if(orientation == "head-to-tail") {
#     cl <- cluster_dimer(d = d, a=50, b=20, c=20, dihedral = 0)
#     res <- spectrum_dispersion(cl, Incidence = pi/2, Axes="x", material=material)
#     return(res[res$polarisation == "p", ])
#   }
#   if(orientation == "side-by-side"){
#     cl <- cluster_dimer(d = d, a=20, b=20, c=50, dihedral = 0)
#     res <- spectrum_dispersion(cl, Incidence = pi/2, Axes="x", material=material)
#     return(res[res$polarisation == "s", ])
#   }
# }

UnitQuaternion(0, 0, 1, 0) # rotation π/2 about y
UnitQuaternion(√2/2, √2/2, 0, 0)
rotation_angle(UnitQuaternion(√2/2, √2/2, 0, 0))
rotation_axis(UnitQuaternion(√2/2, √2/2, 0, 0))
RotZYZ(π/2,0,0)
RotMatrix(UnitQuaternion(√2/2, √2/2, 0, 0))

RotY(π/2)


function model(; d=80, orientation="head-to-tail")
    if orientation == "head-to-tail"
        # dimer axis along y, particle axis along y
        cl = cluster_dimer(d, 20., 40, 20, 0)
        ## incidence: along z (no rotation)
        Incidence = [RotZ(0.)]
        res = spectrum_dispersion(cl, mat, Incidence)
        d = dispersion_df(res, mat.wavelengths)
        return(filter(:polarisation => ==("p"), d))
    elseif orientation == "side-by-side"
        # dimer axis along y, particle axis along z
        cl = cluster_dimer(d, 20., 20, 40, 0)
        ## incidence: along x so rotate z around y axis by π/2
        Incidence = [RotY(π/2)]
        res = spectrum_dispersion(cl, mat, Incidence)
        d = dispersion_df(res, mat.wavelengths)
        return(filter(:polarisation => ==("s"), d))
    end
    
    error("unknown orientation given")
end


wavelength = collect(400:2:800.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)
# cl = cluster_dimer(80, 20., 40, 20, 0)
# Incidence = [RotZ(0.),RotZ(π/2.),RotZ(0.),RotZ(π/2.),RotZ(0.)]
# res = spectrum_dispersion(cl, mat, Incidence)
# d = dispersion_df(res, mat.wavelengths)
# print(d)

params = expand_grid(d=range(80, 200, step=20), orientation=("head-to-tail", "side-by-side"))

all = pmap_df(params, p -> model(; p...))


s = spectrum_dispersion(cl0, mat, [RotZ(0.)])
single = dispersion_df(s, mat.wavelengths)

using AlgebraOfGraphics, CairoMakie

# set_aog_theme!()


xy = data(all) * mapping(:wavelength, :value, color=:d => nonnumeric, col=:orientation, row=:crosstype, linestyle=:crosstype => nonnumeric)
layer = visual(Lines)
draw(layer * xy, facet=(; linkyaxes=:none))


x = repeat(LinRange(-5, 5, 100), inner=2)
df = DataFrame(x=x, y=x .^ 2, z=x .^ 3, col=repeat(1:2, outer=100))

df |>
@vlplot(facet = {column = {field = :col, typ = :nominal}}) +
(
    @vlplot(x = :x) +
    @vlplot(:line, y = :y) +
    @vlplot(mark = {:line, color = :black}, strokeDash = {value = [2, 2]}, y = :z)
)

single |>
@vlplot(facet = {column = {field = :variable, typ = :nominal},
    row = {field = :type, typ = :nominal}}) +
(
    @vlplot() +
    @vlplot(mark = {:line}, data = all, row = :type, column = :variable,
        encoding = {x = "wavelength:q", y = "value:q", color = "d:n", strokeDash = "ϕ:n"}
    ) +
    @vlplot(mark = {:line, color = :black}, data = single,
        strokeDash = {value = [2, 2]}, row = :type, column = :variable,
        encoding = {x = "wavelength:q", y = "value:q"}
    )
)

@vlplot(mark = {:line}, data = all, row = :type, column = :variable,
    encoding = {x = "wavelength:q", y = "value:q", color = "d:n", strokeDash = "ϕ:n"}
)


using VegaDatasets
dataset("weather.csv") |>
@vlplot(repeat = {column = [:temp_max, :precipitation, :wind]}) +
(
    @vlplot() +
    @vlplot(
        :line,
        y = {field = {repeat = :column}, aggregate = :mean, typ = :quantitative},
        x = "month(date):o",
        detail = "year(date):t",
        color = :location,
        opacity = {value = 0.2}
    ) +
    @vlplot(
        :line,
        y = {field = {repeat = :column}, aggregate = :mean, typ = :quantitative},
        x = "month(date):o",
        color = :location
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
grp = mapping(color=:Cyl => nonnumeric);
scat = visual(Scatter)
pipeline = cols * scat
data(mpg) * pipeline |> draw