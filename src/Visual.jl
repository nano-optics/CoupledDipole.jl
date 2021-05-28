function visualise(cl)

        s = []

        for i ∈ 1:length(cl.positions)
                x = cl.positions[i][1]
                y = cl.positions[i][2]
                z = cl.positions[i][3]

                a = cl.sizes[i][1]
                b = cl.sizes[i][2]
                c = cl.sizes[i][3]

                # old Euler formulation
                # α = pi/2+cl.angles[i][1]
                # β = cl.angles[i][2]
                # γ = cl.angles[i][3]
                # particle$i.rotation.z = $α;
                # particle$i.rotation.x = $β;
                # particle$i.rotation.y = $γ;
                # particle$i.rotation.order = 'ZXY';

                # now from quaternion
                # var qm = new THREE.Quaternion();
                # particle$i.useQuaternion = true;
                # particle$i.applyQuaternion();
                # particle$i.quaternion.normalize();
                # or .setRotationFromAxisAngle( axis : Vector3, angle : Float )
                # note: now active rotations so need inverse
                q = Rotations.params(inv(cl.rotations[i]))

                # testing fixed angles
                # alpha=pi/2
                # q = Rotations.params(UnitQuaternion(cos(alpha/2),  0,0,sin(alpha/2)))
                # q = Rotations.params(UnitQuaternion(cos(alpha/2),  0,0,sin(alpha/2)))
                q1 = q[1]; q2 = q[2]; q3 = q[3]; q4 = q[4]
                # note threejs orders them differently (last)
                push!(s,
              "var particle$i = createSphere(radius, segments);
                   particle$i.position.x = $x;
                   particle$i.position.y = $y;
                   particle$i.position.z = $z;
                   particle$i.scale.set($a,$b,$c);
                   particle$i.useQuaternion = true;
                   var qm = new THREE.Quaternion($q2,$q3,$q4,$q1);
                   particle$i.applyQuaternion(qm);
                   particle$i.quaternion.normalize();
                   scene.add(particle$i);\n",
                )
# The X axis is red. The Y axis is green. The Z axis is blue
        end

        s


end
