# Conventions


"""
    <div id="webgl"></div>
    <script>
      (function () {
     
     	var webglEl = document.getElementById('webgl');
     
     	if (!Detector.webgl) {
     		Detector.addGetWebGLMessage(webglEl);
     		return;
     	}
     
     	var width  = window.innerWidth,
     		height = window.innerHeight;
     
     	// sphere params
     	var radius   = 0.2,
     		segments = 36,
     		rotation = 6;  
     
     	var scene = new THREE.Scene();
     
     	var camera = new THREE.PerspectiveCamera(20, width / height, 0.01, 1000);
     	
     	camera.position.z = 10;
     
     	var renderer = new THREE.WebGLRenderer();
     	renderer.setSize(width, height);
     
     	scene.add(new THREE.AmbientLight(0xffffff, 1));
     
     	var xDistance = 0.5;
         var yDistance = 0.5;
         //initial offset so does not start in middle.
         var xOffset = 0;
         var yOffset = 0.2;
         
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
             		
     
     	var controls = new THREE.TrackballControls(camera);
     
     	webglEl.appendChild(renderer.domElement);
     
     	render();
     
     	function render() {
     		controls.update();
     		requestAnimationFrame(render);
     		renderer.render(scene, camera);
     	}
     
     	function createSphere(radius, segments) {
     		return new THREE.Mesh(
     			new THREE.SphereGeometry(radius, segments, segments),
     		);
     	}
     
     
     }());
    </script>
"""

