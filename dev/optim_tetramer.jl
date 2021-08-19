# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using Rotations
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using ForwardDiff



"""
    tetramer of identical spheres, with 3 fixed lengths (d)
    3 free parameters x1..x3
    1st at 0,0,0
    2nd at d,0,0
    3rd at d*cos(x1), d*cos(x1), 0
    4th at d*cos(x2)*sin(x3), d*sin(x2)*sin(x3), d*cos(x3)

"""
function cluster_tethered(x; r = 5, d = 15, material = "Au", type = "particle")
    sizes = [SVector(r, r, r) for _ ∈ 1:4] # identical spheres
    angles = [zero(UnitQuaternion) for _ ∈ 1:4] # angles irrelevant for spheres
    positions = [
        SVector(0.0, 0.0, 0.0),
        SVector(d, 0.0, 0.0),
        SVector(d*cos(x[1]), d*sin(x[1]), 0.0),
        SVector(d*cos(x[2])*sin(x[3]), d*sin(x[2])*sin(x[3]), d*cos(x[3])),
    ]

    Cluster(positions, angles, sizes, [material for _ ∈ 1:4], type)
end

"""
    tetramer of identical spheres, 6 params
    1st at 0,0,0
    2nd at x1,0,0
    3rd at x2, y2, 0
    4th at x3,y3,z3

"""
function cluster_tetramer(x; r = 5, material = "Au", type = "particle")
    sizes = [SVector(r, r, r) for _ ∈ 1:4] # identical particles
    angles = [zero(UnitQuaternion) for _ ∈ 1:4] # angles irrelevant
    positions = [
        SVector(0.0, 0.0, 0.0),
        SVector(x[1], 0.0, 0.0),
        SVector(x[2], x[3], 0.0),
        SVector(x[4], x[5], x[6]),
    ]

    Cluster(positions, angles, sizes, [material for _ ∈ 1:4], type)
end

"""
    model(p)
From parameters p, create the desired physical model
Here:  circular dichroism spectrum with 4 spherical particles
free parameters are 6 coordinates
"""
function model_tetramer(x)

    wavelength = collect(400:5:800.0)
    media = Dict([("Au", epsilon_Au), ("medium", _ -> 1.33)])
    mat = Material(wavelength, media)
    cl = cluster_tetramer(x)

    oa = spectrum_oa(cl, mat)

    return oa.dichroism.extinction ./ oa.average.extinction # g-factor
end

p0 = [50.0, 40.0, 40.0, 30.0, 30.0, 30.0]

p0 = [50.0, 0.0, 50.0, 0.0, -40, 50]

cl = cluster_tetramer(p0)
pl = visualise_makie(cl)
using GLMakie
GLMakie.activate!()
display(pl)

dich = model_tetramer(p0)
 sum(abs.(dich))

 function model_tethered(x)

     wavelength = collect(400:5:800.0)
     media = Dict([("Au", epsilon_Au), ("medium", _ -> 1.33)])
     mat = Material(wavelength, media)
     cl = cluster_tethered(x)

     oa = spectrum_oa(cl, mat)

     return oa.dichroism.extinction ./ oa.average.extinction # g-factor
 end


 p0 = [pi/4, pi/3, pi/5]

 using CairoMakie
 CairoMakie.activate!()
 cl = cluster_tethered(p0)
 pl = visualise_makie(cl)
 display(pl)

 dich = model_tethered(p0)
  sum(abs.(dich))

## testing the model works as expected

using Base.Iterators


d1 = map(x -> DataFrame(wavelength = collect(400:5:800.0),
                    value = model([50.0, 0.0, 50.0, 0.0, 50.0, x]), position = x, group="z"), 10:10:60)

d1 = map(x -> DataFrame(wavelength = collect(400:5:800.0),
                        value = model([50.0, 0.0, 50.0, 0.0, x, 50]), position = x, group="x"), -60:10:60)

# d2 = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
#                         value = model([0.1,π/180 * θ,π/4]), angle=θ, group="β"), 0:15:90)
#
# d3 = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
#                         value = model([π/180 * θ,0.1,0.5]), angle=θ, group="γ"), 0:15:90)

# m = vcat(d1...,d2...,d3...)
m = vcat(d1...)


@vlplot(data=m,
    mark = "line",
    row = "group",
    encoding = {x = "wavelength:q", y = "value:q",color="position:n"}
)


combine(groupby(m, [:angle, :group]), :value => x -> sum(abs.(x)))

# function to minimise
function objective_dichroism(x; radius=50)

    # test for collisions by checking all pairwise distances
    positions = hcat(
        [0.0, 0.0, 0.0],
        [0.0, x[1], 0.0],
        [x[2], x[3], 0.0],
        [x[4], x[5], x[6]])

    distances = Distances.pairwise(Euclidean(), positions, dims=2)

     if(any(vec(distances) .< radius))
         return -Inf
     end

    dichroism = model(x)
    return -sum(abs.(dichroism))
end
#
# objective_dichroism([0,π/4,0])
#
# gradFD = ForwardDiff.gradient(objective_dichroism, [0.0,0,0])
#
# # first, 1D problem
#
# x0 = [0.1]
# function cost_1d(x)
#     -sum(abs.(model([0,0,x])))
# end
#
# using Optim
#
# pf = optimize(cost_1d, 0, π/2, Brent())
# # pf.minimizer
# # 0.7857593236287779
# # π/4
# # 0.7853981633974483
#
#
# # using ReverseDiff
# # gradRD = ReverseDiff.gradient(objective_dichroism, [0.0,0,0])
# # # super slow and returns NaN
#
#
# x0 = [0.01,0.01,0.7]
# pf = optimize(objective_dichroism, x0, LBFGS(), Optim.Options(iterations = 100);
#          autodiff = :forward)
#
# pf.minimizer
# 3-element Vector{Float64}:
#  -6.8631915491679e-17
#   0.01
#   0.7857593353326985


# using Zygote

# gradZY = Zygote.gradient(objective_dichroism, [0.0,0,0])
# errors?
