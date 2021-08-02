# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using ForwardDiff


function cluster_two(α_1, α_2, β; d=100, a=20, b=20, c=30, material = "Au", type="particle")
    sizes = [SVector(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector(0.0, y, 0.0) for y in (-d/2, d/2)]
    angles = [SVector(α_1, β, 0),
              SVector(α_2, 0, 0)]
    Cluster(positions, angles, sizes, material, type)
end


"""
    model(p)
From parameters p, create the desired physical model
Here:  circular dichroism spectrum with 2 spheroidal particles at fixed positions
free parameters are 3 angles
"""
function model(p)

    wavelength = collect(400:5:800.0)
    media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
    mat = Material(wavelength, media)
    cl = cluster_two(p[1], p[2], p[3])
    oa = spectrum_oa(cl, mat)

    return oa.dichroism.extinction
end

# dich = model([0,π/4,0])
# sum(abs.(dich))


## testing the model works as expected

using Base.Iterators
# collect(product([1,2,3], [4,5,6]))

# p1 = [(0,0,α) for α ∈ 0:15:90]
# p2 = [(0,β,π/4) for β ∈ 0:15:90]
# p3 = [(γ,0.1,0.5) for γ ∈ 0:15:90]

d1 = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
                    value = model([0,0,π/180 * θ]), angle=θ, group="α"), 0:15:90)

d2 = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
                        value = model([0.1,π/180 * θ,π/4]), angle=θ, group="β"), 0:15:90)

d3 = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
                        value = model([π/180 * θ,0.1,0.5]), angle=θ, group="γ"), 0:15:90)

m = vcat(d1...,d2...,d3...)


@vlplot(data=m,
    mark = "line",
    row = "group",
    encoding = {x = "wavelength:q", y = "value:q",color="angle:n"}
)


combine(groupby(m, [:angle, :group]), :value => x -> sum(abs.(x)))


# function to minimise
function objective_dichroism(p)
    dich = model(p)
    return -sum(abs.(dich))
end

objective_dichroism([0,π/4,0])

gradFD = ForwardDiff.gradient(objective_dichroism, [0.0,0,0])

# first, 1D problem

x0 = [0.1]
function cost_1d(x)
    -sum(abs.(model([0,0,x])))
end

using Optim

pf = optimize(cost_1d, 0, π/2, Brent())
# pf.minimizer
# 0.7857593236287779
# π/4
# 0.7853981633974483


# using ReverseDiff
# gradRD = ReverseDiff.gradient(objective_dichroism, [0.0,0,0])
# # super slow and returns NaN


x0 = [0.01,0.01,0.7]
pf = optimize(objective_dichroism, x0, LBFGS(), Optim.Options(iterations = 100);
         autodiff = :forward)

pf.minimizer
# 3-element Vector{Float64}:
#  -6.8631915491679e-17
#   0.01
#   0.7857593353326985


# using Zygote

# gradZY = Zygote.gradient(objective_dichroism, [0.0,0,0])
# errors?
