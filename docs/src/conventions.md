# Conventions


```@raw html
    <div id="webgl"></div>
    <script>
      (function () {

	var webglEl = document.getElementById('webgl');

	if (!Detector.webgl) {
		Detector.addGetWebGLMessage(webglEl);
		return;
	}

	var width  = 400,
		height = 400;

	// sphere params
	var radius   = 1,
		segments = 36,
		rotation = 6;  

	var scene = new THREE.Scene();
scene.background = new THREE.Color( 0xeeeeee );


	var camera = new THREE.PerspectiveCamera(20, width / height, 0.01, 1000);
	
	camera.position.z = 10;

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
         
        var i=0,j=0;
        var sphere1 = createSphere(radius, segments);
        sphere1.position.y = - yOffset;
        sphere1.position.x = - xOffset;
        sphere1.scale.set( 0.5, 0.5, 0.8 );
        scene.add(sphere1);
        
        var sphere2 = createSphere(radius, segments);
        sphere2.position.y = yOffset;
        sphere2.position.x = xOffset;
        sphere2.scale.set( 0.5, 0.5, 0.8 );
        
       // sphere2.rotation.x = 1.57;
        sphere2.rotation.y = 0.7853982;
        scene.add(sphere2);
        		
//var axesHelper = new THREE.AxesHelper( 5 );
//scene.add( axesHelper );

				scene.add( new THREE.AxesHelper( 1 ) );
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


}());
    </script>
```

