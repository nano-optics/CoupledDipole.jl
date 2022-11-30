# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using Rotations

## this example looks at a helical strand of Au nanorods in water
## contrasting angular dispersion of CD against OA

## materials
wavelengths = collect(450:2:750.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

# cl <- cluster_helix(5, R0=20, pitch=30,
#                     delta=pi/2, delta0=0, right=TRUE,
#                     a=15/2, b=15/2, c=15,
#                     angles="helix")

# nn <- 7
# Incidence <- rep(seq(0,pi/2,length=nn), 3)


## dimer geometry
# cluster_helix(N, a, b, c, R, Λ, δ, δ_0=0, handedness="left",
#     material="Au", type="particle")
# d, a, b, c, ϕ=0.0
# cl1 = cluster_helix(5, 15.0 / 2, 15.0 / 2, 15.0, 20, 30, π / 2, 0)
cl1 = cluster_helix(5, 15.0 / 2, 15.0 / 2, 15.0, 20, 25, π / 2, 0)
cl0 = cluster_single(15.0 / 2, 15.0 / 2, 15.0)


# visualise_threejs(cl1)


## incidence: along z, along x, along y
angles = range(0, π / 2, 7)
Incidence = QuatRotation.([RotX.(angles); RotY.(angles); RotZ.(angles)])

disp1 = spectrum_dispersion(cl0, mat, Incidence, polarisations="circular")
disp2 = spectrum_dispersion(cl1, mat, Incidence, polarisations="circular")

dd1 = dispersion_df(disp1, mat.wavelengths, format="wide")
dd2 = dispersion_df(disp2, mat.wavelengths, format="wide")

rotations = DataFrame(:angle => repeat(angles * 180 / π, outer=3),
      :axis => repeat(["x", "y", "z"], inner=length(angles)),
      :angle_id => eachindex(Incidence))


# tmp = stack(leftjoin(dd2, rotations, on=:angle_id), [:pol1, :pol2], variable_name=:type)
# tmp2 = stack(leftjoin(dd1, rotations, on=:angle_id), [:pol1, :pol2], variable_name=:type)

# plt = data(filter(:crosstype => ==("absorption"), tmp)) * 
#     mapping(:wavelength, :value, color=:angle => nonnumeric, col=:axis,
#       row=:type) * visual(Lines)

# plt2 = data(filter(:crosstype => ==("absorption"), tmp2)) * 
# mapping(:wavelength, :value, color=:angle => nonnumeric, col=:axis,
#   row=:type) * visual(Lines, linestyle=:dash)
# draw(plt + plt2, facet=(; linkyaxes=:none),
#       palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))


oa0 = spectrum_oa(cl0, mat)
oa1 = spectrum_oa(cl1, mat)

d0 = oa_df(oa0, mat.wavelengths)
d3 = oa_df(oa1, mat.wavelengths)

using DataFramesMeta
dich = @chain dd2 begin
      @rsubset :crosstype == "extinction"
      @transform :dichroism = :pol1 .- :pol2
      @transform :average = 0.5 .* (:pol1 .+ :pol2)
end

rotations = DataFrame(:angle => repeat(angles * 180 / π, outer=3),
      :axis => repeat(["x", "y", "z"], inner=length(angles)),
      :angle_id => eachindex(Incidence))

toplot = leftjoin(dich, rotations, on=:angle_id)
d4 = stack(toplot, [:dichroism, :average], variable_name=:type)

using ColorSchemes
set_aog_theme!()
d2 = data(filter(:crosstype => ==("extinction"), d3))
d1 = data(d4)
m1 = d1 * mapping(:wavelength, :value, color=:angle => nonnumeric, col=:axis,
      row=:type)
m2 = d2 * mapping(:wavelength, :value, row=:type)
layer1 = m1 * visual(Lines)
layer2 = m2 * visual(Lines, linestyle=:dash)
draw(layer1 + layer2, facet=(; linkyaxes=:none),
      palettes=(; color=cgrad(ColorSchemes.phase.colors, 12, categorical=true)))


# d5 = [insertcols!(d3, :cluster => "helix", :hand => "left")
#       insertcols!(d0, :cluster => "single", :hand => "_")]


# d5 |> @vlplot(
#       width = 120,
#       height = 100,
#       mark = {:line},
#       row = "type", column = "crosstype",
#       resolve = {scale = {y = "independent"}},
#       encoding = {x = "wavelength:q", y = "value:q",
#             strokeDash = "cluster:n", color = "hand:n"}
# )
