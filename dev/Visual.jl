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

cl =cluster_helix(10, 20, 20, 40, 50, 300)


p = visualise(cl)


io = open("myfile.txt", "w");
write(io, broadcast(*, p...));
close(io);

function euler_to_axisangle(e)
        M = euler_passive(e...)
        θ = acos((M[1, 1] + M[2, 2] + M[3, 3] - 1) / 2)+1e-5
        e1 = (M[3, 2] - M[2, 3]) / (2sin(θ))
        e2 = (M[1, 3] - M[3, 1]) / (2sin(θ))
        e3 = (M[2, 1] - M[1, 2]) / (2sin(θ))

        θ .* Vec3f0(e1, e2, e3+1e-5)
end

# euler_to_axisangle([0,pi/2,0])

cl = cluster_dimer(10, 1, 2, 3, pi/4)

meshscatter(
    Point3f0.(cl.positions),
    marker =sphere,
    markersize = 1.0.*Vec3f0.(cl.sizes),
    #markersize = [Vec3f0(10,10,20) for i in 1:length(cl.positions)],
    rotation = euler_to_axisangle.(cl.angles),
    color=:red, limits = Rect([0,0,0],[1,1,1])
)







using Makie
meshscatter(
    rand(Point3f0, 10),
    marker =sphere,
    markersize = [Vec3f0(0.1,0.1,0.2) for i in 1:10],
    rotation = [Vec3f0(0,  0,pi/2) for i in 1:10],
    color=:red
)





using Makie
using GLMakie
p1 = Point3(0.0, 0.0, 0.0)
p2 = Point3(0.0, 1.0, 0.0)
p3 = Point3(0.0, 1.0, 1.0)
p4 = Point3(0.0, 0.0, 1.0)

sc = mesh([p1, p2, p3], color = :blue, shading = false)
