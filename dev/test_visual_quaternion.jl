# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite

cl = cluster_helix(10, 10, 10, 20, 50, 200)


p = visualise(cl)

io = open("myfile.txt", "w");
write(io, broadcast(*, p...));
close(io);
