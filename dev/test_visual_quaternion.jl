# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite
using Rotations

# cl = cluster_helix(10, 10, 10, 20, 50, 200)

# cl = cluster_dimer(50, 10,10,50,π/4, 20*pi/180,10*pi/180)

function cluster_cone(N, a, b, c, ϕ = 0.0, α = 0.0)
    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    positions = [SVector(0.0, 0.0, 0.0) for i in 1:N]
    q2 = UnitQuaternion(cos(α/2), sin(α/2), 0, 0) # rotation α about x
    # rotate particle 1 by q1 only (stays in yz plane)
    # rotate particle 2 by q2, then q3 but in original frame so order swapped
    rotations = [q2 * UnitQuaternion(cos(ϕ/2), 0, sin(ϕ/2), 0)  for ϕ in LinRange(0,360,N+1)*pi/180]
    Cluster(positions, rotations[2:end], sizes, "Au", "particle")
end


cl = cluster_cone(18, 5,5,50,π/4, 60*pi/180)

p = visualise(cl)

io = open("myfile.txt", "w");
write(io, broadcast(*, p...));
close(io);
