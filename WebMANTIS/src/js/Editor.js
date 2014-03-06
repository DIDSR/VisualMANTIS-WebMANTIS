/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        webMANTIS v1.0		     //
//			     //	     (optical photons transport visualization)       //
//			     //							     //
//			     //////////////////////////////////////////////////////////
//
// 
//
//
// ****Disclaimer****
//  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
//  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
//  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
//  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
//  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
//  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
//  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
//  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
//  derived from it, and any modified versions bear some notice that they have been modified. 
//
//
//
//	Filename:	Editor.js
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Signals are callbacks when the user perform a specific function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
var Editor = function () {

	var SIGNALS = signals;

	this.signals = {

		// notifications
		rendererChanged: new SIGNALS.Signal(),
		sceneGraphChanged: new SIGNALS.Signal(),
		windowResize: new SIGNALS.Signal(),
		
		//WebMANTIS related callbacks
		executeStarted: new SIGNALS.Signal(),
		renderScene: new SIGNALS.Signal(),
		nextHistory: new SIGNALS.Signal(),
		saveInputArguments: new SIGNALS.Signal()
	};


	this.scene = new THREE.Scene();
	this.sceneHelpers = new THREE.Scene();

	this.object = {};

	this.selected = null;
	this.helpers = {};

};

Editor.prototype = {

	saveInputArguments: function()
	{
		this.signals.saveInputArguments.dispatch();
	},
	
	//signal for when user presses Execute
	executeStart: function()
	{
		this.signals.executeStarted.dispatch();
	},
	
	setScene: function ( scene ) {

		this.scene.name = scene.name;
		this.scene.userData = JSON.parse( JSON.stringify( scene.userData ) );

		while ( scene.children.length > 0 ) {

			this.addObject( scene.children[ 0 ] );

		}

	},

	// adds the spheres/optical photons to the scene
	addSpheres: function()
	{
		//sets previous visible to false
		for(var i = 0; i < spheres.length; i++)
		{
			for(var j=0;j < spheres[i].length; j++)
			{
				spheres[i][j].visible = false;
			}
		}
		for(var i = 0; i < ceilings.length; i++)
		{
			for(var j=0;j < ceilings[i].length; j++)
			{
				ceilings[i][j].visible = false;
			}
		}
		for(var i = 0; i < floors.length; i++)
		{
			for(var j=0;j < floors[i].length; j++)
			{
				floors[i][j].visible = false;
			}
		}
		for(var i = 0; i < lines.length; i++)
		{
			for(var j=0;j < lines[i].length; j++)
			{
				lines[i][j].visible = false;
			}
		}
		
		//sets up the different colors for the spheres/optical photons
		var colorArray = new Array();
		
		var red = new THREE.MeshBasicMaterial( { color: 0xFF0000, transparent: true, opacity: 1.0 } );
		var green = new THREE.MeshBasicMaterial( { color: 0x00FF00, transparent: true, opacity: 1.0 } );
		var blue = new THREE.MeshBasicMaterial( { color: 0x0000FF, transparent: true, opacity: 1.0 } );
		var yellow = new THREE.MeshBasicMaterial( { color: 0xFFFF00, transparent: true, opacity: 1.0} );
		var aqua = new THREE.MeshBasicMaterial( { color: 0x00FFFF, transparent: true, opacity: 1.0} );
		var pink = new THREE.MeshBasicMaterial( { color: 0xFF00FF, transparent: true, opacity: 1.0} );
		var grey = new THREE.MeshBasicMaterial( { color: 0xC0C0C0, transparent: true, opacity: 1.0} );
		var orange = new THREE.MeshBasicMaterial( { color: 0xFF6600, transparent: true, opacity: 1.0} );
		var purple = new THREE.MeshBasicMaterial( { color: 0x9900FF, transparent: true, opacity: 1.0} );
		
		colorArray.push(red);
		colorArray.push(green);
		colorArray.push(blue);
		colorArray.push(yellow);
		colorArray.push(aqua);
		colorArray.push(pink);
		colorArray.push(grey);
		colorArray.push(aqua);
		colorArray.push(orange);
		colorArray.push(purple);
		
		//sets up the different colors for the trajectories/lines
		var lineColorArray = new Array();
		
		var lred = new THREE.LineBasicMaterial( { color: 0xFF0000 } );
		var lgreen = new THREE.LineBasicMaterial( { color: 0x00FF00 } );
		var lblue = new THREE.LineBasicMaterial( { color: 0x0000FF } );
		var lyellow = new THREE.LineBasicMaterial( { color: 0xFFFF00 } );
		var laqua = new THREE.LineBasicMaterial( { color: 0x00FFFF } );
		var lpink = new THREE.LineBasicMaterial( { color: 0xFF00FF } );
		var lgrey = new THREE.LineBasicMaterial( { color: 0xC0C0C0 } );
		var lorange = new THREE.LineBasicMaterial( { color: 0xFF6600 } );
		var lpurple = new THREE.LineBasicMaterial( { color: 0x9900FF } );
		
		lineColorArray.push(lred);
		lineColorArray.push(lgreen);
		lineColorArray.push(lblue);
		lineColorArray.push(lyellow);
		lineColorArray.push(laqua);
		lineColorArray.push(lpink);
		lineColorArray.push(lgrey);
		lineColorArray.push(laqua);
		lineColorArray.push(lorange);
		lineColorArray.push(lpurple);
		
		//sets up the floor and ceiling boxes
		var pgeometry = new THREE.PlaneGeometry( 20, 20, 1, 1);
		var pmaterialFloor = new THREE.MeshBasicMaterial( { color: 0xFFFFFF, transparent: true, opacity: 1.0 } );
		var pmaterialCeil = new THREE.MeshBasicMaterial( { color: 0x000000, transparent: true, opacity: 1.0 } );
		
		
		var k = 0;
		var sphereColor = colorArray[k];
		var lineColor = lineColorArray[k];
		
		var tempArray = new Array();
		var ceilTempArray = new Array();
		var floorTempArray = new Array();
		
		var lineTempArray = new Array();
		var sphere, oldSphere, newSphere;
		
		var lineGeometry = new THREE.Geometry();
		
		var pointCounter = 0;
					
		//reads in the JSON data returned by the PHP script
		$.each(globalJSONData, function(key, val) 
		{
			//gets the x, y, z coordinates
			var zz = parseFloat(val.z) + (parseFloat(detector2) / 2.0);//75.0;
			var xx = parseFloat(val.x - ( parseFloat(xdim2) / 2.0)); //(909.0 / 2));
			var yy = parseFloat(val.y - ( parseFloat(ydim2) / 2.0)); //(909.0 / 2));
			
			// adds the spheres and lines connecting spheres to the scene
			if(zz >= 0 && zz <= parseFloat(detector2))//150.0)
			{
				if(zz == 0)
				{
					var pmeshFloor = new THREE.Mesh( pgeometry, pmaterialFloor );
					pmeshFloor.rotation.x = - Math.PI/2;
					pmeshFloor.position.set(xx, zz, yy);
					floorTempArray.push(pmeshFloor);
				}
				else if(zz == parseFloat(detector2))//150.0)
				{
					var pmeshCeil = new THREE.Mesh( pgeometry, pmaterialCeil );
					pmeshCeil.rotation.x = - Math.PI/2;
					pmeshCeil.position.set(xx, zz, yy);
					ceilTempArray.push(pmeshCeil);
				}
				
				sphere = new THREE.Mesh(new THREE.SphereGeometry(1), sphereColor);
				sphere.position.set(xx, zz, yy);
				tempArray.push(sphere);
				
				//adds the line to connect two spheres together
				if(pointCounter == 0)
				{
					pointCounter += 1;
					oldSphere = sphere;
					newSphere = sphere;
				}
				else
				{
					pointCounter += 1;
					oldSphere = newSphere;
					newSphere = sphere;
					
					lineGeometry = new THREE.Geometry();
					
					lineGeometry.vertices.push(oldSphere.position);
					lineGeometry.vertices.push(newSphere.position);	
					
					var line = new THREE.Line( lineGeometry, lineColor );
					lineTempArray.push(line);				
				}
				
				//changes to a new color when one history ends
				if(parseInt(val.terminate) > 0)
				{
					k = k + 1;
					sphereColor = colorArray[k];
					lineColor = lineColorArray[k];
					pointCounter = 0;
				}
			}
			else if (parseInt(val.terminate) > 0)
			{
				k = k + 1;
				sphereColor = colorArray[k];
				lineColor = lineColorArray[k];
				pointCounter = 0;
			}
		});
		
		//pushes the data into global data structures
		lines.push(lineTempArray);
		spheres.push(tempArray);
		ceilings.push(ceilTempArray);
		floors.push(floorTempArray);
		
		//sets new data visibility to true
		for(var i=0;i<spheres[spheres.length-1].length;i++)
		{
			this.scene.add(spheres[spheres.length-1][i]);
			spheres[spheres.length-1][i].visible = true;
		}
		
		for(var i=0;i<ceilings[ceilings.length-1].length;i++)
		{
			this.scene.add(ceilings[ceilings.length-1][i]);
			ceilings[ceilings.length-1][i].visible = true;
		}
		
		for(var i=0;i<floors[floors.length-1].length;i++)
		{
			this.scene.add(floors[floors.length-1][i]);
			floors[floors.length-1][i].visible = true;
		}
		
		for(var i=0;i<lines[lines.length-1].length;i++)
		{
			this.scene.add(lines[lines.length-1][i]);
			lines[lines.length-1][i].visible = true;
		}
				
		this.signals.renderScene.dispatch();
		this.signals.nextHistory.dispatch();
	},

	contains: function (a, obj)
	{
		for (var i = 0; i < a.length; i++) 
		{
	        if (a[i] === obj) 
	        {
	            return true;
	        }
	    }
	    return false;
	},
	
	// Code to add cylinders of the detectorto the scene
	addCylinder: function ()
	{
		//sets previous cylinders to false
		for(var i = 0; i < objects.length;i++)
		{
			for(var j=0;j<objects[i].length;j++)
			{
				objects[i][j].visible = false;
			}
		}
		
		//sets up the colors
		var colorArray = new Array();
		
		var cylinderGeometry = new THREE.CylinderGeometry( parseFloat(colRad2), parseFloat(colRad2), parseFloat(detector2), 20, 4); //5.1, 5.1, 150, 20, 4 );
		var red = new THREE.MeshBasicMaterial( { color: 0xFF0000, transparent: true, opacity: 0.1 } );
		var green = new THREE.MeshBasicMaterial( { color: 0x00FF00, transparent: true, opacity: 0.1 } );
		var blue = new THREE.MeshBasicMaterial( { color: 0x0000FF, transparent: true, opacity: 0.1} );
		var yellow = new THREE.MeshBasicMaterial( { color: 0xFFFF00, transparent: true, opacity: 0.1} );
		var aqua = new THREE.MeshBasicMaterial( { color: 0x00FFFF, transparent: true, opacity: 0.1} );
		var pink = new THREE.MeshBasicMaterial( { color: 0xFF00FF, transparent: true, opacity: 0.1} );
		var grey = new THREE.MeshBasicMaterial( { color: 0xC0C0C0, transparent: true, opacity: 0.1} );
		var orange = new THREE.MeshBasicMaterial( { color: 0xFF6600, transparent: true, opacity: 0.1} );
		var purple = new THREE.MeshBasicMaterial( { color: 0x9900FF, transparent: true, opacity: 0.1} );
		
		colorArray.push(red);
		colorArray.push(green);
		colorArray.push(blue);
		colorArray.push(yellow);
		colorArray.push(aqua);
		colorArray.push(pink);
		colorArray.push(grey);
		colorArray.push(aqua);
		colorArray.push(orange);
		colorArray.push(purple);
		
		var tempArray = new Array();
		var k = 0;
		var cylinderMaterial = colorArray[k];
				
		//for each cylinder in the JSON returned by the PHP script
		$.each(globalJSONData, function(key, val) 
		{			
			//creates a new cylinder and sets its location
			var cylinder = new THREE.Mesh( cylinderGeometry, cylinderMaterial );
			cylinder.position.set(parseFloat(val.Xc - (parseFloat(xdim2) / 2)), (parseFloat(detector2) / 2.0), parseFloat(val.Yc - (parseFloat(ydim2) / 2)));//cylinder.position.set(parseFloat(val.Xc - (909.0 / 2)), 75.0, parseFloat(val.Yc - (909.0 / 2)));
			
			var bool = false;
			for (var i = 0; i < tempArray.length; i++) 
			{
				if (tempArray[i].position.x == cylinder.position.x && 
						tempArray[i].position.y == cylinder.position.y) 
				{
					bool = true;
				}
			}
	    
			if(bool == false)
			{
				tempArray.push(cylinder);
			}
			
			// move to the next color if the current history ends
			if(parseInt(val.terminate) > 0)
			{
				k = k + 1;
				cylinderMaterial = colorArray[k];
				
				if(k == colorArray.length)
				{
					alert("Too many colors used, resetting to 0");
					k = 0;
				}
			}
		});
		
		objects.push(tempArray);
		for(var i=0;i<objects[objects.length-1].length;i++)
		{
			this.scene.add(objects[objects.length-1][i]);
			objects[objects.length-1][i].visible = true;
		} 
		
		this.signals.renderScene.dispatch();
	}
}
