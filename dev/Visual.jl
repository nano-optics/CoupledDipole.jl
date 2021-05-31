# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using Rotations

cl = cluster_dimer(10, 1, 2, 3, pi/4)

cl = cluster_helix(8, 20, 20, 40, 100, 300, π/4, 0, "right")


function as_quat(v)
 Quaternionf0(v[1],v[2],v[3],v[4])
  Vec4f0(v...)
end

function as_rotation(q)
  v = rotation_axis(q)
  qrotation(Vec3f0(v[1],v[2],v[3]), rotation_angle(q))
end



function Makie_rotation(q)
        q = inv(q) # passify
        v = Rotations.rotation_axis(q)
        θ = Rotations.rotation_angle(q)
        qrotation(Vec3f0(v...), θ)
end


function visualise_makie(cl; colour=:gold)

        meshscatter(
                Point3f0.(cl.positions),
                markersize = Vec3f0.(cl.sizes),
                rotations = Makie_rotation.(cl.rotations),
                color = colour,
        )

end

using GLMakie
GLMakie.activate!()

visualise_makie(cl, colour = :gold)



meshscatter(
    Point3f0.(cl.positions),
    markersize = Vec3f0.(cl.sizes),
    rotations = as_rotation.(inv.(cl.rotations)),
    color=:gold
)



function axisxangle(q)
Vec3f0(rotation_angle((q)) .* rotation_axis((q)))
end



using GLMakie

xs = cos.(1:0.5:20)
ys = sin.(1:0.5:20)
zs = LinRange(0, 3, length(xs))

meshscatter(xs, ys, zs, markersize = [0.1,0.2,0.5], color = zs)




using Makie
using GLMakie
p1 = Point3(0.0, 0.0, 0.0)
p2 = Point3(0.0, 1.0, 0.0)
p3 = Point3(0.0, 1.0, 1.0)
p4 = Point3(0.0, 0.0, 1.0)

sc = mesh([p1, p2, p3], color = :blue, shading = false)

function visualise(cl)

        s = []

        for i ∈ 1:length(cl.positions)
                x = cl.positions[i][1]
                y = cl.positions[i][2]
                z = cl.positions[i][3]

                a = cl.sizes[i][1]
                b = cl.sizes[i][2]
                c = cl.sizes[i][3]

                α = pi/2+cl.angles[i][1]
                β = cl.angles[i][2]
                γ = cl.angles[i][3]
                push!(s,
                        "var sphere$i = createSphere(radius, segments);
                   sphere$i.position.x = $x;
                   sphere$i.position.y = $y;
                   sphere$i.position.z = $z;
                   sphere$i.scale.set($a,$b,$c);
                   sphere$i.rotation.z = $α;
                   sphere$i.rotation.x = $β;
                   sphere$i.rotation.y = $γ;
                   sphere$i.rotation.order = 'ZXY';
                   scene.add(sphere$i);\n",
                )
        end

        s


end

cl = cluster_helix(10, 10, 10, 20, 50, 200)


p = visualise(cl)

io = open("myfile.txt", "w");
write(io, broadcast(*, p...));
close(io);



using WebIO
Asset("https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js")



function euler_to_axisangle(e)
        M = euler_passive(e...)
        θ = acos((M[1, 1] + M[2, 2] + M[3, 3] - 1) / 2)+1e-5
        e1 = (M[3, 2] - M[2, 3]) / (2sin(θ))
        e2 = (M[1, 3] - M[3, 1]) / (2sin(θ))
        e3 = (M[2, 1] - M[1, 2]) / (2sin(θ))

        θ .* Vec3f0(e1, e2, e3+1e-5)
end

# euler_to_axisangle([0,pi/2,0])
