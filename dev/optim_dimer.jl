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
    dimer of identical spheroids, with fixed gap (d)
    3 free parameters x1..x3
    1st at 0,0,0
    2nd at d,0,0
    3rd at d*cos(x1), d*cos(x1), 0
    4th at d*cos(x2)*sin(x3), d*sin(x2)*sin(x3), d*cos(x3)

"""
function cluster_dimer3(
    x,
    a = 10.0,
    b = 10.0,
    c = 20.0,
    d = 80.0,
    material = "Au",
    type = "particle",
)
    cluster_dimer(
        d,
        a,
        b,
        c,
        x[1],
         x[2],
         x[3],
         material,
         type
    )
end



# cluster_dimer(80, 10, 10, 20, 1,2,3, "Au","part")

 function model_dimer(x; wavelength = collect(400:5:800.0))

     media = Dict([("Au", epsilon_Au), ("medium", _ -> 1.33)])
     mat = Material(wavelength, media)
     cl = cluster_dimer3(x)

     oa = spectrum_oa(cl, mat)

          # return oa.dichroism.extinction ./ oa.average.extinction # g-factor
     return oa.dichroism.extinction
 end


 p0 = [0.7, -0.2, 0.1]

cluster_dimer3(p0,1)

 using CairoMakie
 CairoMakie.activate!()
 cl = cluster_dimer3(p0)
 pl = visualise_makie(cl)
 display(pl)

 dich = cluster_dimer3(p0)

  dich = model_dimer(p0)

  sum(abs.(dich))

## testing the model works as expected

using Base.Iterators


d1 = map(x -> DataFrame(wavelength = collect(400:2:800.0),
                    value = model_dimer([x*pi/180, 0.1,0.2], wavelength = collect(400:2:800.0)), dihedral = x, group="z"), 0:15:90)

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
    encoding = {x = "wavelength:q", y = "value:q",color="dihedral:n"}
)


combine(groupby(m, [:angle, :group]), :value => x -> sum(abs.(x)))

# function to minimise
function objective_dichroism(x; wavelength = collect(550:1:560.0))

    dichroism = model_dimer(x, wavelength=wavelength)
    return -sum(abs.(dichroism))
    # return -max(abs.(dichroism)...)
end


objective_dichroism([π/4,-π/4,π/2])
objective_dichroism([0,π/4,0])
#
gradFD = ForwardDiff.gradient(objective_dichroism, [π/4,-π/4,π/2])

using Optim

x0 = [0.7,0.02,0.01]
pf = optimize(objective_dichroism, x0, LBFGS(), Optim.Options(iterations = 100);
         autodiff = :forward)
#
pf.minimizer
# 3-element Vector{Float64}:




clf = cluster_dimer3(pf.minimizer)
pl = visualise_makie(clf)
display(pl)


mi = DataFrame(wavelength = collect(400:2:800.0), type = "initial",
                    value = model_dimer(p0, wavelength = collect(400:2:800.0)))

mf = DataFrame(wavelength = collect(400:2:800.0), type = "final",
                    value = model_dimer(pf.minimizer, wavelength = collect(400:2:800.0)))

m = vcat(mi, mf)


@vlplot(data=m,
    mark = "line",
    encoding = {x = "wavelength:q", y = "value:q", color = "type:n"}
)
