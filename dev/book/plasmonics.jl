# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using Rotations
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using DataFramesMeta
using VegaLite
using AlgebraOfGraphics, CairoMakie
using ColorSchemes
set_aog_theme!()

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

params = expand_grid(a0=[10, 30, 50], ar=[1, 1.2, 1.5, 2], material=["Au", "Si"])
all = pmap_df(params, p -> model(; p...))

all
map = mapping(:wavelength, :value, row=:crosstype,
    col=:a0 => nonnumeric => "a0", linestyle=:material, color=:ar => nonnumeric => "ar")
l1 = data(@rsubset(all, :type == "average")) * map * visual(Lines)
draw(l1, facet=(; linkyaxes=:none))

