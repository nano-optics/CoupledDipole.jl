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

# UnitQuaternion(0, 0, 1, 0) # rotation π/2 about y
# UnitQuaternion(√2/2, √2/2, 0, 0)
# rotation_angle(UnitQuaternion(√2/2, √2/2, 0, 0))
# rotation_axis(UnitQuaternion(√2/2, √2/2, 0, 0))
# RotZYZ(π/2,0,0)
# RotMatrix(UnitQuaternion(√2/2, √2/2, 0, 0))

# RotY(π/2)

# initial version, but now prefer low level cluster
# function model(; d=80, orientation="head-to-tail")
#     if orientation == "head-to-tail"
#         # dimer axis along y, particle axis along y
#         cl = cluster_dimer(d, 20., 40, 20, 0)
#         ## incidence: along z (no rotation)
#         Incidence = [RotZ(0.)]
#         res = spectrum_dispersion(cl, mat, Incidence)
#         d = dispersion_df(res, mat.wavelengths)
#         return(filter(:polarisation => ==("p"), d))
#     elseif orientation == "side-by-side"
#         # dimer axis along y, particle axis along z
#         cl = cluster_dimer(d, 20., 20, 40, 0)
#         ## incidence: along x so rotate z around y axis by π/2
#         Incidence = [RotY(π/2)]
#         res = spectrum_dispersion(cl, mat, Incidence)
#         d = dispersion_df(res, mat.wavelengths)
#         return(filter(:polarisation => ==("s"), d))
#     end

#     error("unknown orientation given")
# end

function model(; d=100, orientation="head-to-tail")
    if orientation == "head-to-tail"
        # dimer axis along y
        positions = [SVector(0.0, y, 0.0) for y in (-d / 2.0, d / 2.0)]
    elseif orientation == "side-by-side"
        # dimer axis along x
        positions = [SVector(x, 0.0, 0.0) for x in (-d / 2.0, d / 2.0)]
    end

    sizes = [SVector(20.0, 50.0, 20.0), SVector(20.0, 50.0, 20.0)]
    rotations = repeat([UnitQuaternion(1.0, 0.0, 0.0, 0.0)], 2)
    materials = repeat(["Au"], 2)
    cl = Cluster(positions, rotations, sizes, materials, "particle")

    Incidence = [RotZ(0.0)] ## incidence: along z (no rotation)
    res = spectrum_dispersion(cl, mat, Incidence)
    d = dispersion_df(res, mat.wavelengths)
end

## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)


params = expand_grid(d=range(100, 500, step=50), orientation=("head-to-tail", "side-by-side"))
all = pmap_df(params, p -> model(; p...))


## reference geometry
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

cl0 = cluster_single(20.0, 50.0, 20.0)
# cl0 = cluster_single(50.0, 20.0, 20.0)
# s = spectrum_dispersion(cl0, mat, [UnitQuaternion(RotZ(0.0))])
s = spectrum_dispersion(cl0, mat, [SVector(0.0, 0.0, 0.0)])
single = dispersion_df(s, mat.wavelengths)


xy = data(single) * mapping(:wavelength, :value, row=:crosstype, col=:polarisation)
layer = visual(Lines)
draw(layer * xy, facet=(; linkyaxes=:none))


# d = [insertcols(all, :cluster => "dimer"),
#     insertcols(single, :cluster => "single", :d => missing, :orientation => missing)]

# @vlplot(data = d,
#     width = 400,
#     height = 300,
#     mark = {:line},
#     row = "crosstype",
#     resolve = {scale = {y = "independent"}},
#     encoding = {x = "wavelength:q", y = "value:q", color = "d:n", strokeDash = "cluster:n"}
# )

using AlgebraOfGraphics, CairoMakie

# set_aog_theme!()


xy = data(filter(:polarisation => ==("p"), all)) * mapping(:wavelength, :value, color=:d => nonnumeric, col=:orientation, row=:crosstype)
layer = visual(Lines)
draw(layer * xy, facet=(; linkyaxes=:none))






# s = spectrum_oa(cl0, mat)
# single = oa_df(s, mat.wavelengths)


# xy = data(single) * mapping(:wavelength, :value, color=:type, row=:crosstype,col=:type)
# layer = visual(Lines)
# draw(layer * xy, facet=(; linkyaxes=:none))

using Arrow
Arrow.write("test.arrow", single)
Arrow.write("testc.arrow", single, compress=:lz4)


y = Arrow.Table("test.arrow") |> DataFrame
x = Arrow.Table("testr.arrow") |> DataFrame
eltype.(eachcol(x))

# using CSV

# CSV.write("test.csv", single)



