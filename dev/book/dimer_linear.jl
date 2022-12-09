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

## this example looks at 2 Au nanorods in water
## contrasting head-t-t and side-b-s configurations
## extinction spectra for varying distances
## fixed orientation

## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

function model(; d=100, orientation="head-to-tail")
    if orientation == "head-to-tail"
        # dimer axis along y
        positions = [SVector(0.0, y, 0.0) for y in (-d / 2.0, d / 2.0)]
    elseif orientation == "side-by-side"
        # dimer axis along x
        positions = [SVector(x, 0.0, 0.0) for x in (-d / 2.0, d / 2.0)]
    end

    sizes = [SVector(20.0, 50.0, 20.0), SVector(20.0, 50.0, 20.0)]
    rotations = repeat([QuatRotation(1.0, 0.0, 0.0, 0.0)], 2)
    materials = repeat(["Au"], 2)
    cl = Cluster(positions, rotations, sizes, materials, "particle")

    Incidence = [RotZ(0.0)] ## incidence: along z (no rotation)
    res = spectrum_dispersion(cl, mat, Incidence)
    d = dispersion_df(res, mat.wavelengths)
end


params = expand_grid(d=range(100, 500, step=50), orientation=("head-to-tail", "side-by-side"))
all = pmap_df(params, p -> model(; p...))


## reference geometry
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

sizes = [SVector(0.0, 0.0, 0.0)]
positions = [SVector(0.0, 0.0, 0.0)]
# input parameters are Euler angles
rotations = [QuatRotation(RotZYZ(0.0, 0.0, 0.0))]
Cluster(positions, rotations, sizes, ["Au"], "particle")

cl0 = cluster_single(20.0, 50.0, 20.0)
s = spectrum_dispersion(cl0, mat, [QuatRotation(RotZ(0.0))])
single = dispersion_df(s, mat.wavelengths)


# xy = data(single) * mapping(:wavelength, :value, row=:crosstype, col=:polarisation)
# layer = visual(Lines)
# draw(layer * xy, facet=(; linkyaxes=:none))

d1 = data(@rsubset(all, :polarisation == "pol2"))
d2 = data(@rsubset(single, :polarisation == "pol2"))

m1 = d1 * mapping(:wavelength, :value, color=:d => nonnumeric, col=:orientation, row=:crosstype)
m2 = d2 * mapping(:wavelength, :value, row=:crosstype)
layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
draw(layer1 + layer2, facet=(; linkyaxes=:none), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))
# https://docs.juliaplots.org/latest/generated/colorschemes/
# cgrad(ColorSchemes.viridis.colors, 12, categorical=true)
# cgrad([:purple, :green], 12, categorical=true)
# ColorSchemes.viridis.colors
# cgrad(ColorSchemes.viridis.colors, 12, categorical=true)

# using Arrow
# Arrow.write("test.arrow", single)

# y = Arrow.Table("test.arrow") |> DataFrame
# x = Arrow.Table("testr.arrow") |> DataFrame

