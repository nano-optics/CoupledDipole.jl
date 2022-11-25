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
    q2 = QuatRotation(cos(α/2), sin(α/2), 0, 0) # rotation α about x
    # rotate particle 1 by q1 only (stays in yz plane)
    # rotate particle 2 by q2, then q3 but in original frame so order swapped
    rotations = [q2 * QuatRotation(cos(ϕ/2), 0, sin(ϕ/2), 0)  for ϕ in LinRange(0,360,N+1)*pi/180]
    Cluster(positions, rotations[2:end], sizes, ["Au" for ii in 1:N], ["particle" for ii in 1:N])
end


cl = cluster_cone(18, 5,5,50,π/4, 60*pi/180)

 cl = cluster_helix(13, 10, 10, 20, 50, 200)
p = visualise(cl)

io = open("mycluster.html", "w");
write(io, p);
close(io);




using StaticArrays
using WebIO

struct Stuff{T}
   positions::Vector{SVector{3,T}}
   sizes::Vector{SVector{3,T}}
   angles::Vector{SVector{3,T}}
end


cl = Stuff([@SVector rand(3) for _ in 1:4],
[@SVector rand(3) for _ in 1:4],
[@SVector rand(3) for _ in 1:4])

function WebIO.render(cl::Stuff)
    return dom"div"(
        # dom"p"(header),
        dom"p"(cl.positions),
        # dom"p"(footer)    ,
    )
end

cl
WebIO.render(cl)


using WebIO
w = Scope(imports=["//cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"])
w(dom"div.container"(dom"p"("this is a test")))





dom"script"(js"""test""")


using WebIO, JSExpr

w = Scope(imports=["//cdnjs.cloudflare.com/ajax/libs/p5.js/0.5.11/p5.js"])
onmount(w, @js function (p5)
    function sketch(s)
        s.setup = () -> s.createCanvas(640, 480)

        s.draw = function ()
          if s.mouseIsPressed
            s.fill(0)
          else
            s.fill(255)
          end
          s.ellipse(s.mouseX, s.mouseY, 80, 80)
        end
    end
    @new p5(sketch, this.dom.querySelector(".container"))
end)
w(dom"div.container"())
