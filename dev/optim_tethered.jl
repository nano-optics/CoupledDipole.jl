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

 function model_tethered(x; d=60, r=20, wavelength = collect(400:5:800.0))

     media = Dict([("Au", epsilon_Au), ("medium", _ -> 1.33)])
     mat = Material(wavelength, media)
     cl = cluster_tethered(x, d=d, r=r)

     oa = spectrum_oa(cl, mat)

     return oa.dichroism.extinction ./ oa.average.extinction # g-factor
 end


 p0 = [pi/4, pi/3, 1.96]

 using CairoMakie
 CairoMakie.activate!()
 cl = cluster_tethered(p0)
 pl = visualise_makie(cl)
 display(pl)

 dich = model_tethered(p0)

  dich = model_tethered(p0, d=20)

  sum(abs.(dich))

## testing the model works as expected

using Base.Iterators


d1 = map(x -> DataFrame(wavelength = collect(400:5:800.0),
                    value = model_tethered([pi/4, pi/3, x]), position = x, group="z"), -pi:pi/8:pi)

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

using Distances

function pair_distances(positions) # upper triangular distances
    A = Distances.pairwise(Distances.Euclidean(), positions, dims=2)
   return [A[i, j] for j in 2:size(A,1) for i in 1:j-1]
end

# function to minimise
function objective_dichroism(x; r = 20, d = 80, wavelength = collect(550:1:560.0))

    # test for collisions by checking all pairwise distances
    positions = hcat(
        [0.0, 0.0, 0.0],
        [d, 0.0, 0.0],
        [d * cos(x[1]), d * sin(x[1]), 0.0],
        [d * cos(x[2]) * sin(x[3]), d * sin(x[2]) * sin(x[3]), d * cos(x[3])],
    )

    if (any(pair_distances(positions) .< 2.5 * r))
        return +Inf
    end

    dichroism = model_tethered(x, d = d, r = r, wavelength=wavelength)
    # return -sum(abs.(dichroism))
    return -max(abs.(dichroism)...)
end

#
cl = cluster_tethered([π/4,-π/4,π/2], d=80,r=20)
pl = visualise_makie(cl)
display(pl)

d = 80
x = [π/4,-π/4,π/2]
positions = hcat(
    [0.0, 0.0, 0.0],
    [d, 0.0, 0.0],
    [d * cos(x[1]), d * sin(x[1]), 0.0],
    [d * cos(x[2]) * sin(x[3]), d * sin(x[2]) * sin(x[3]), d * cos(x[3])],
)

pair_distances(positions)

objective_dichroism([π/4,-π/4,π/2])
objective_dichroism([0,π/4,0])
#
gradFD = ForwardDiff.gradient(objective_dichroism, [π/4,-π/4,π/2])

using Optim

x0 = [π/4,-π/4,π/2]
pf = optimize(objective_dichroism, x0, LBFGS(), Optim.Options(iterations = 100);
         autodiff = :forward)
#
pf.minimizer
# 3-element Vector{Float64}:




clf = cluster_tethered(pf.minimizer, d=80,r=20)
pl = visualise_makie(clf)
display(pl)

x = pf.minimizer
positions = hcat(
    [0.0, 0.0, 0.0],
    [d, 0.0, 0.0],
    [d * cos(x[1]), d * sin(x[1]), 0.0],
    [d * cos(x[2]) * sin(x[3]), d * sin(x[2]) * sin(x[3]), d * cos(x[3])],
)
pair_distances(positions)


m = DataFrame(wavelength = collect(400:5:800.0),
                    value = model_tethered(pf.minimizer, d=60,r=20))

@vlplot(data=m,
    mark = "line",
    encoding = {x = "wavelength:q", y = "value:q"}
)
