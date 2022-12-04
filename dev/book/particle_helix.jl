# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using DataFramesMeta
using VegaLite
using AlgebraOfGraphics
using Makie
using Rotations
using ColorSchemes
set_aog_theme!()

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

disp0 = spectrum_dispersion(cl0, mat, Incidence, polarisations="circular")
disp1 = spectrum_dispersion(cl1, mat, Incidence, polarisations="circular")

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
lay2 = dat2 * mapping(:wavelength, :value, color=:angle => nonnumeric, col=:axis,
             row=:type) * visual(Lines)
draw(lay1 + lay2, facet=(; linkyaxes=:none),
      palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))


