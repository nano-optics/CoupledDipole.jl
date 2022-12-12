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

params = expand_grid(radius=range(10.0, 50, step=10), gap=range(5.0, 20, step=5))
all = pmap_df(params, p -> model(p...))

using Arrow
terms_d = Arrow.Table("./dev/book/terms_data.arrow") |> DataFrame
# terms_d[:,:crosstype] .= a 
terms_d[:,:cross] = replace.(terms_d[:,:crosstype], "Abs" => "absorption", "Sca" => "scattering", "Ext" => "extinction")


d1 = data(@rsubset(all, :type == "average",:crosstype == "extinction"))
d2 = data(@rsubset(terms_d, :dipole == "1",:cross == "extinction"))
d3 = data(@rsubset(terms_d, :dipole == "Inf",:cross == "extinction"))

# m1 = d1 * mapping(:wavelength, :value, color=:gap => nonnumeric => "gap /nm", col=:crosstype, row=:radius => nonnumeric)

# m2 = d2 * mapping(:wavelength, :average => x -> x ./ 4, color=:gap => nonnumeric => "gap /nm", 
# col=:cross, row=:radius => nonnumeric)

m1 = d1 * mapping(:wavelength, :value, col=:gap => nonnumeric => "gap /nm", row=:radius => nonnumeric)

m2 = d2 * mapping(:wavelength, :average => x -> x ./ 4, col=:gap => nonnumeric => "gap /nm", 
row=:radius => nonnumeric)

m3 = d3 * mapping(:wavelength, :average => x -> x ./ 4, col=:gap => nonnumeric => "gap /nm", 
row=:radius => nonnumeric)



layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
layer3 = m3 * visual(Lines, linestyle=:dot)
fg=draw(layer1 + layer2 + layer3, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))



save("figure.pdf", fg, px_per_unit=3)

## look at dichroism


d1 = data(@rsubset(all, :type == "dichroism",:crosstype == "extinction"))
d2 = data(@rsubset(terms_d, :dipole == "1",:cross == "extinction"))
d3 = data(@rsubset(terms_d, :dipole == "Inf",:cross == "extinction"))


m1 = d1 * mapping(:wavelength, :value, col=:gap => nonnumeric => "gap /nm", row=:radius => nonnumeric)

m2 = d2 * mapping(:wavelength, :average => x -> x ./ 4, col=:gap => nonnumeric => "gap /nm", 
row=:radius => nonnumeric)

m3 = d3 * mapping(:wavelength, :average => x -> x ./ 4, col=:gap => nonnumeric => "gap /nm", 
row=:radius => nonnumeric)



layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
layer3 = m3 * visual(Lines, linestyle=:dot)
fg=draw(layer1 + layer2 + layer3, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))



save("figure2.pdf", fg, px_per_unit=3)

