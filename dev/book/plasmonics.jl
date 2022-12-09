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
using AlgebraOfGraphics, Makie, CairoMakie
home = homedir()
const font_folder = "$home/Library/Fonts/"
firasans(weight) = joinpath(font_folder, "FiraSans-$(weight).ttf")
cmu(weight) = joinpath(font_folder, "cmun$(weight).ttf")
# set_aog_theme!(fonts=[cmu("rm"), cmu("rm")])

gill(weight) = joinpath(font_folder, "GillSansNova-$(weight).otf")
set_aog_theme!(fonts=[gill("Book"), gill("Light")])

# update_theme!(                  # Tweaks the current theme
#     fonts=[firasans("Regular"), firasans("Light")],
#     fontsize=12
# )
# 

## this example illustrates plasmon resonances in metal nanoparticles

## materials

function model(; a0, ar, material="Au", wavelength=collect(450:2:850.0))
    # options for different materials
    tab_fun = Dict([("Au", epsilon_Au), ("Ag", epsilon_Ag), ("Si", epsilon_Si)])
    epsilon_fun = tab_fun[material]
    media = Dict([(material, epsilon_fun), ("medium", x -> 1.33)])
    mat = Material(wavelength, media)

    a, c = spheroid_ar(a0, ar)
    cl = cluster_single(a, a, c, 0.0, 0.0, 0.0, material)

    oa = spectrum_oa(cl, mat)
    d = oa_df(oa, mat.wavelengths)
    # d[!, :material] .= material
    # d[!, :ar] .= ar
    # d[!, :a0] .= a0
    return d
end


test = model(; a0=30.0, ar=1.5)

params = expand_grid(a0=[10, 20, 30, 40], ar=[1, 1.5, 2], material=["Au"])
gold = pmap_df(params, p -> model(; p...))

map = mapping(:wavelength, :value, col=:crosstype,
    row=:a0 => nonnumeric => "a0", color=:ar => nonnumeric => "ar",
    linestyle=:material)
l1 = data(@rsubset(gold, :type == "average")) * map * visual(Lines)
draw(l1, facet=(; linkyaxes=:none))

# params = expand_grid(a0=[10, 20, 30, 40], ar=[1, 1.5, 2], material=["Ag"])

# silver = pmap_df(params, p -> model(; p...,(wavelength=collect(350:2:850.0))))

# map = mapping(:wavelength, :value, col=:crosstype,
#     row=:a0 => nonnumeric => "a0", color=:ar => nonnumeric => "ar",
#     linestyle=:material)
# l1 = data(@rsubset(silver, :type == "average")) * map * visual(Lines)
# draw(l1, facet=(; linkyaxes=:none))

## non metal
params2 = expand_grid(a0=[10, 40], ar=[1, 1.5, 2], material=["Si"])
silicon = pmap_df(params2, p -> model(; p...))

map = mapping(:wavelength, :value, col=:crosstype,
    row=:a0 => nonnumeric => "a0", color=:ar => nonnumeric => "ar")
l1 = data(@rsubset(silicon, :type == "average")) * map * visual(Lines)
draw(l1, facet=(; linkyaxes=:none), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"))

