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
    Cluster(positions, rotations[2:end], sizes, ["Au" for ii in 1:N], ["particle" for ii in 1:N])
end


cl = cluster_cone(18, 5,5,50,π/4, 60*pi/180)

 cl = cluster_helix(13, 10, 10, 20, 50, 200)
p = visualise(cl)

io = open("mycluster.html", "w");
write(io, p);
close(io);


header = js"""
<div id='webgl'></div>

	<script src='https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js'></script>
	<script src='https://unpkg.com/three@0.85.0/examples/js/controls/TrackballControls.js'></script>
	<script src='https://unpkg.com/three@0.85.0/examples/js/Detector.js'></script>

	<script>
// Created by Bjorn Sandvik - thematicmapping.org
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
         var yOffset = 1;
     """


footer = js"""
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
	// material = new THREE.MeshBasicMaterial( {color: 0xffff00} );
//	 const material = new THREE.MeshLambertMaterial({
//  color: 0xdaa520,
//  emissive: 0x111111,
//});
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


}());"""


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
        dom"p"(header),
        dom"p"(cl.positions),
        dom"p"(footer)    ,
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
