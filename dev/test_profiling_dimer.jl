# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite

## materials
wavelength = collect(450:2:850.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, media)

## dimer geometry
cl1 = cluster_dimer(80, 20, 20, 40, π / 4)
cl2 = cluster_single(20, 20, 40)

## incidence: along z, along x, along y
Incidence = [SVector(0, 0, 0), SVector(0, π / 2, 0), SVector(π / 2, π / 2, 0)]

@time spectrum_dispersion(cl1, mat, Incidence)
@time spectrum_oa(cl1, mat)

using Profile
@profile spectrum_dispersion(cl1, mat, Incidence)
# Juno.@profiler spectrum_dispersion(cl1, mat, Incidence)

disp1 = spectrum_dispersion(cl1, mat, Incidence)
disp2 = spectrum_dispersion(cl2, mat, Incidence)

d1 = dispersion_df(disp1, mat.wavelength)
d2 = dispersion_df(disp2, mat.wavelength)

d = [insertcols!(d1, :cluster => "dimer")
    insertcols!(d2, :cluster => "single")]

@vlplot(data = d,
    width = 400,
    height = 300,
    mark = {:line},
    row = "crosstype",
    resolve = {scale = {y = "independent"}},
    encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash = "cluster:n"}
)



oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl2, mat)


# Juno.@enter spectrum_oa(cl1, mat)



# aa = DataFrame(:wavelength => mat.wavelength, :value => oa1.average.absorption)
# @vlplot(data=aa,
# width= 400,
# height =  300,
#     mark = {:line},
#     encoding = {x = "wavelength:q", y = "value:q"}
# )
#




d3 = oa_df(oa1, mat.wavelength)
d4 = oa_df(oa2, mat.wavelength)

d5 = [insertcols!(d3, :cluster => "dimer")
    insertcols!(d4, :cluster => "single")]

# filter(row -> row[:type] == "average" , d5)
d5 |> @vlplot(
    width = 400,
    height = 300,
    mark = {:line},
    row = "type",
    resolve = {scale = {y = "independent"}},
    encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash = "cluster:n"}
)
