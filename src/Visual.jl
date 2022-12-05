"""
     visualise_threejs(cl; colour=:gold)

Standalone interactive 3D visual of a cluster

Examples

```
cl = cluster_helix(8, 20, 20, 40, 100, 300, π/4, 0, "right")
visualise_threejs(cl)
# open cluster.html in browser
```

"""
function visualise_threejs(cl; out="cluster.html")

        header = "<!doctype html>
<html>
<head>
	<meta charset='utf-8'>
	<meta name='viewport' content='width=device-width, initial-scale=1.0'>
	<title>Threejs cluster visualisation</title>
	<style>
		body { margin: 0; overflow: hidden; background-color: #fff; }
		.tm  { position: absolute; top: 10px; right: 10px; }
		.webgl-error { font: 15px/30px monospace; text-align: center; color: #fff; margin: 50px; }
		.webgl-error a { color: #fff; }
    </style>
</head>
<body>
        <div id='webgl'></div>

        	<script src='https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js'></script>
        	<script src='https://unpkg.com/three@0.85.0/examples/js/controls/TrackballControls.js'></script>
        	<script src='https://unpkg.com/three@0.85.0/examples/js/Detector.js'></script>

        	<script>

        (function () {

        	var webglEl = document.getElementById('webgl');

        	if (!Detector.webgl) {
        		Detector.addGetWebGLMessage(webglEl);
        		return;
        	}

        	var width  = 600,
        		height = 600;

        	// sphere params
        	var radius   = 1,
        		segments = 36,
        		rotation = 6;

        	var scene = new THREE.Scene();
        scene.background = new THREE.Color( 0xeeeeee );


        	var camera = new THREE.PerspectiveCamera(100, width / height, 0.01, 1000);

        	camera.position.z = 200;

        	var renderer = new THREE.WebGLRenderer();
        	renderer.setSize(width, height);

        const frontSpot = new THREE.SpotLight(0xeeeece);
        frontSpot.position.set(1000, 1000, 1000);
        scene.add(frontSpot);

        const frontSpot2 = new THREE.SpotLight(0xddddce);
        frontSpot2.position.set(-500, -500, -500);
        scene.add(frontSpot2);


        	scene.add(new THREE.AmbientLight(0xffffff, 1));

             	var xDistance = 2;
               var yDistance = 2;
                 //initial offset so does not start in middle.
                 var xOffset = 0;
                 var yOffset = 1;\n
             "


        footer = "
        	scene.add( new THREE.AxesHelper( 100 ) );
        	var controls = new THREE.TrackballControls(camera);

        	webglEl.appendChild(renderer.domElement);

        	render();

        	function render() {
        		controls.update();
        		requestAnimationFrame(render);
        		renderer.render(scene, camera);
        	}



        	function createSphere(radius, segments) {

        const material = new THREE.MeshStandardMaterial({
          color: 0xfcc742,
          emissive: 0x111111,
          specular: 0xffffff,
          metalness: 0.8,
          roughness: 0.6,
        });
        		return new THREE.Mesh(
        			new THREE.SphereGeometry(radius, segments, segments),
        			material
        		);
        	}


        }());
</script>
</body>
</html>\n"

        s = [header]

        for i ∈ 1:length(cl.positions)
                x = cl.positions[i][1]
                y = cl.positions[i][2]
                z = cl.positions[i][3]

                a = cl.sizes[i][1]
                b = cl.sizes[i][2]
                c = cl.sizes[i][3]

                # note: active rotations so need inverse
                q = Rotations.params(inv(cl.rotations[i]))

                # testing fixed angles
                # alpha=pi/2
                # q = Rotations.params(QuatRotation(cos(alpha/2),  0,0,sin(alpha/2)))
                q1 = q[1]
                q2 = q[2]
                q3 = q[3]
                q4 = q[4]
                # note threejs orders them differently (last)
                push!(
                        s,
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

        push!(s, footer)

        io = open(out, "w")
        write(io, broadcast(.*, s...))
        close(io)

end



function Makie_rotation(q)
        q = inv(q) # passify
        v = Rotations.rotation_axis(q)
        θ = Rotations.rotation_angle(q)
        qrotation(Vec3f0(v...), θ)
end


"""
     visualise_makie(cl; colour=:gold)

Static 3D visual of a cluster

Examples

```
cl = cluster_helix(8, 20, 20, 40, 100, 300, π/4, 0, "right")
visualise_makie(cl, colour = :silver)
```

"""
function visualise_makie(cl; colour=:gold)

        meshscatter(
                Point3f0.(cl.positions),
                markersize=Vec3f0.(cl.sizes),
                rotations=Makie_rotation.(cl.rotations),
                color=colour,
        )

end
