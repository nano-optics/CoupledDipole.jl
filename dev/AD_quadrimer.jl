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
# using BenchmarkTools


function cluster_two(α, β_1, β_2; d=100, a=20, b=20, c=30, material = "Au", type="particle")
    sizes = [SVector(a, b, c) for ii in 1:2] # identical particles
    positions = [SVector(0.0, y, 0.0) for y in (-d/2, d/2)]
    angles = [SVector(α, β_1, 0),
              SVector(0, β_2, 0)]
    Cluster(positions, angles, sizes, material, type)
end


"""
    physics(p)
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


## testing the model works as expected

using Iterators
# collect(product([1,2,3], [4,5,6]))

d = map(θ -> DataFrame(wavelength = collect(400:5:800.0),
                    cd = model([0,θ*pi/180,0]), theta=θ), 0:15:90)

m = vcat(d...)

@vlplot(data=m,
width= 400,
height =  300,
    mark = {:line},
    encoding = {x = "wavelength:q", y = "cd:q", color = "theta:n"}
)



dich = model([0,π/4,0])
# sum(abs.(dich))


# function to minimise
function objective_dichroism(p)
    dich = model(p)
    return -sum(abs.(dich))
end

objective_dichroism([0,π/4,0])

gradFD = ForwardDiff.gradient(objective_dichroism, [0.0,0,0])

# using ReverseDiff
#
# gradRD = ReverseDiff.gradient(objective_dichroism, [0.0,0,0])
# # super slow and returns NaN


using Optim


x0 = [0.0,0.1,0.2]
optimize(objective_dichroism, x0, LBFGS(); autodiff = :forward)


# using Zygote

# gradZY = Zygote.gradient(objective_dichroism, [0.0,0,0])
# errors?
