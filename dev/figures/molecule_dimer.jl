# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
# using VegaLite
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

## materials
wavelength = collect(400:1:700.0)
media = Dict([("Rhodamine", alpha_lorentzmolecule), ("medium", x -> 1.33)])
mat = Material(wavelength, media)


## dimer geometry
cl0 = cluster_single(0, 0, 1, 0, 0, 0, "Rhodamine", "point")
cl1 = cluster_dimer(0.8, 0, 0, 1, 0, 0, 0, "Rhodamine", "point")
cl2 = cluster_dimer(0.8, 0, 0, 1, π / 4, 0, 0, "Rhodamine", "point")

## incidence: along z, along x, along y
Incidence = [RotZ(0.0), RotY(π / 2), RotX(π / 2)]

disp1 = spectrum_dispersion(cl0, mat, Incidence)
disp2 = spectrum_dispersion(cl1, mat, Incidence)

d1 = incidence_labels(dispersion_df(disp1, mat.wavelengths), Incidence, ['z', 'y', 'x'])
d2 = incidence_labels(dispersion_df(disp2, mat.wavelengths), Incidence, ['z', 'y', 'x'])


dd1 = data(@rsubset(d1, :polarisation == "pol1"))
dd2 = data(@rsubset(d2, :polarisation == "pol1"))

m = mapping(:wavelength, :value, color=:axis_label => nonnumeric,
    row=:crosstype, col=:axis_label)
l1 = visual(Lines)
l2 = visual(Lines, linestyle=:dash)
draw(dd1 * m * l1 + dd2 * m * l2, facet=(; linkyaxes=:rowwise), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))




# oa1 = spectrum_oa(cl0, mat)
# oa2 = spectrum_oa(cl2, mat)


# d3 = oa_df(oa1, mat.wavelengths)
# d4 = oa_df(oa2, mat.wavelengths)

# d5 = [insertcols!(d4, :cluster => "dimer")
#     insertcols!(d3, :cluster => "single")]

# d5 |> @vlplot(
#     width = 400,
#     height = 300,
#     mark = {:line},
#     row = "type",
#     resolve = {scale = {y = "independent"}},
#     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash = "cluster:n"}
# )
