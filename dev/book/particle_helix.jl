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

## this example looks at a helical strand of Au nanorods in water
## contrasting angular dispersion of CD against OA

## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

# cluster_helix(N, a, b, c, R, Λ, δ, δ_0=0, handedness="left",
#     material="Au", type="particle")
cl1 = cluster_helix(5, 15.0 / 2, 15.0 / 2, 15.0, 20, 25, π / 2, 0)
cl0 = cluster_single(15.0 / 2, 15.0 / 2, 15.0)

# visualise_threejs(cl1)

## incidence: along z, along x, along y
angles = range(0, π / 2, 7)
Incidence = QuatRotation.([RotX.(angles); RotY.(angles); RotZ.(angles)])

disp0 = spectrum_dispersion(cl0, mat, Incidence, polarisation="circular")
disp1 = spectrum_dispersion(cl1, mat, Incidence, polarisation="circular")

disp0 = dispersion_df(disp0, mat.wavelengths, format="wide")
disp1 = dispersion_df(disp1, mat.wavelengths, format="wide")

rotations = DataFrame(:angle => repeat(angles * 180 / π, outer=3),
      :axis => repeat(["x", "y", "z"], inner=length(angles)),
      :angle_id => eachindex(Incidence))

oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)

d0 = oa_df(oa0, mat.wavelengths)
d1 = oa_df(oa1, mat.wavelengths)

dd1 = @chain disp1 begin
      @rsubset :crosstype == "extinction"
      @transform :dichroism = :pol1 .- :pol2
      @transform :average = 0.5 .* (:pol1 .+ :pol2)
end

rotations = DataFrame(:angle => repeat(angles * 180 / π, outer=3),
      :axis => repeat(["x", "y", "z"], inner=length(angles)),
      :angle_id => eachindex(Incidence))

tmp = leftjoin(dd1, rotations, on=:angle_id)
d2 = stack(tmp, [:dichroism, :average], variable_name=:type)

dat1 = data(@rsubset(d1, :crosstype == "extinction"))
dat2 = data(d2)
lay1 = dat1 * mapping(:wavelength, :value, row=:type) * visual(Lines, linestyle=:dash)
lay2 = dat2 * mapping(:wavelength, :value, color=:angle => nonnumeric => "angle /º", col=:axis,
             row=:type) * visual(Lines)
fg = draw(lay1 + lay2, facet=(; linkyaxes=:none), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
      palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))


