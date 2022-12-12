# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFramesMeta
using DataFrames
using Rotations
using LaTeXStrings
using AlgebraOfGraphics, Makie, CairoMakie, ColorSchemes

home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
# set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])

gill(weight) = joinpath(font_folder, "GillSansNova-$(weight).otf")
set_aog_theme!(fonts=[gill("Book"), gill("Light")])


## this example considers a helix of spheres with different separations
## to test the limits of the coupled-dipole approximation


## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

## cluster geometry
function model(radius, gap)
    cl = cluster_quadrimer(radius, gap)
    oa = spectrum_oa(cl, mat)
    return oa_df(oa, mat.wavelengths)
end

## manually test cases

d1 = model(30, 10)
d2 = model(30, 15)
d3 = model(35, 10)

map1 = mapping(:wavelength, :value, row=:type, col=:crosstype)

m1 = map1 * (data(d1) * visual(Lines) +
             data(d2) * visual(Lines, linestyle=:dash) +
             data(d3) * visual(Lines, linestyle=:dot))
fg = draw(m1, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"))

fg

# save("figure.pdf", fg, px_per_unit=3)

## loop

params = expand_grid(radius=range(10, 50, step=10), gap=range(5, 20, step=5))
all = pmap_df(params, p -> model(p...))

d1 = data(@rsubset(all, :type == "average"))

m1 = d1 * mapping(:wavelength, :value, color=:gap => nonnumeric => "gap /nm", col=:crosstype, row=:radius => nonnumeric)

layer1 = m1 * visual(Lines)
# layer2 = m2 * visual(Lines, linestyle=:dash)
draw(layer1, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))

