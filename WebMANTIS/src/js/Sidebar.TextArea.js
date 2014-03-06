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
//	Filename:	Sidebar.TextArea.js
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Loads the side bar of the window
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
Sidebar.TextArea = function ( editor ) 
{
	var signals = editor.signals;

	var container = new UI.Panel();
	container.setPadding( '10px' );
	container.setBorderTop( '1px solid #ccc' );

	//This section sets the number of user accounts
	// has to be manually edited for now
	var userRow = new UI.Panel();
	var wel = new UI.Text( 'User:  ').setColor( '#000' );
	var users = new UI.Select().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).onChange( updateUser );
	users.addOption('Guest');
	//user.addOption('new user');
	userRow.add(wel);
	userRow.add(users);
	container.add(userRow);
	
	// History section
	var topRow = new UI.Panel();
	var percentNow = new UI.Text( 'Current History: 0' ).setColor( '#000' );
	var mypid = new UI.Text('Process ID: 0').setColor('#000');
	
	topRow.add(percentNow);
	topRow.add( new UI.Break());
	topRow.add(mypid);
	
	container.add( topRow );
	container.add (new UI.HorizontalRule());

	// Zip download section
	var zipRow = new UI.Panel();
	
	var saveZip = new UI.Button('Save To Zip').onClick( function ()
	{
		var initialName = currentUser;
		var filename = prompt("Please enter a zip file name:", initialName+".zip");
		
		$.ajax({
			type: "POST",
			url: "tarFile.php",
			data: { folder : currentUser, fname : filename },
			async: false,
			success: function(data)
			{
				var response = jQuery.parseJSON(data); //file is zipped
				alert("Created "+response.filename);
				zipFiles.addOption(response.filename);
			}
		});						
	});
	
	var downloadZip = new UI.Button('Download').onClick( function ()
	{
		var val2 = currentUser+"/"+zipFiles.getValue();
		window.location = val2;
	});
	
	var zipFiles = new UI.Select().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).onChange( updateZip );

	zipRow.add( saveZip );
	zipRow.add( zipFiles );
	zipRow.add( downloadZip );
	container.add(zipRow);
	//container.add( new UI.Break() );
	container.add (new UI.HorizontalRule());
	
	//Output textbox section
	container.add( new UI.Text( 'Output:' ).setColor( '#666' ) );
	container.add( new UI.Break(), new UI.Break() );
	
	var outliner = new UI.TextArea().setWidth( '100%' ).setHeight('140px').setColor( '#444' ).setFontSize( '10px' );
	container.add( outliner );
	container.add( new UI.Break() );
	container.add (new UI.HorizontalRule());

	//Pulse height spectra section
	var phsRow = new UI.Panel();
	var phsLeft = new UI.Button( '⇦' ).onClick( function () 
	{
		phsLeftClick();
	});
	
	var phsRight = new UI.Button( '⇨' ).onClick( function () 
	{
		phsRightClick();
	});
	
	var phsImage = new UI.Canvas('phs');
	
	phsRow.add(phsLeft);
	phsRow.add(phsRight);
	
	container.add(phsImage);
	container.add(phsRow);
	container.add (new UI.HorizontalRule());
	
	//Point response function section
	var prfRow = new UI.Panel();
	var prfLeft = new UI.Button( '⇦' ).onClick( function () 
	{
		prfLeftClick();
	});
	
	var prfRight = new UI.Button( '⇨' ).onClick( function () 
	{
		prfRightClick();
	});
	
	var prfImage = new UI.Canvas('prf');
	
	prfRow.add(prfLeft);
	prfRow.add(prfRight);
	
	container.add(prfImage);
	container.add(prfRow);
	container.add (new UI.HorizontalRule());
	
	container.add( new UI.Text( 'Input Arguments:' ).setColor( '#666' ) );
	container.add( new UI.Break() );
	container.add( new UI.Break() );
	
	//Below are the input arguments
	// X-Ray Histories
	var xrayRow = new UI.Panel();
	var xrayIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('100000').setTitle('Number of X-Ray Histories to be simulated (N)').onChange( update );
	xrayRow.add( new UI.Text( 'X-Ray Histories' ).setWidth( '90px' ).setColor( '#666' ) );
	xrayRow.add( xrayIn );
	container.add( xrayRow );
	
	// Number Photon Histories Visualization
	var photonsRow = new UI.Panel();
	var photonsIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('10').setTitle('Number of photons to be visualized').onChange( update );
	photonsRow.add( new UI.Text( 'Number Photon Histories Visualization' ).setWidth( '90px' ).setColor( '#666' ) );
	photonsRow.add( photonsIn );
	container.add( photonsRow );
	
	// Min Detect
	var mindRow = new UI.Panel();
	var mindIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0').setTitle('Minimum number of optical photons that can be detected').onChange( update );
	mindRow.add( new UI.Text( 'Min Detect' ).setWidth( '90px' ).setColor( '#666' ) );
	mindRow.add( mindIn );
	container.add( mindRow );
	
	// Max Detect
	var maxdRow = new UI.Panel();
	var maxdIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('1400').setTitle('Maximum (N* Yield) number of optical photons that can be detected').onChange( update );
	maxdRow.add( new UI.Text( 'Max Detect' ).setWidth( '90px' ).setColor( '#666' ) );
	maxdRow.add( maxdIn );
	container.add( maxdRow );
	
	// Number of Bins
	var numBinsRow = new UI.Panel();
	var numBinsIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('140').setTitle('Number of bins for storing Pulse Height Spectrum. The maximum value for this paramter is 1000.').onChange( update );
	numBinsRow.add( new UI.Text( 'Number of Bins' ).setWidth( '90px' ).setColor( '#666' ) );
	numBinsRow.add( numBinsIn );
	container.add( numBinsRow );
	
	// X-Dimension
	var xdimRow = new UI.Panel();
	var xdimIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('909').setTitle('X-Dimension of detector (in microns).').onChange( update );
	xdimRow.add( new UI.Text( 'X-Dimension' ).setWidth( '90px' ).setColor( '#666' ) );
	xdimRow.add( xdimIn );
	container.add( xdimRow );
	
	// Y-Dimension
	var ydimRow = new UI.Panel();
	var ydimIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('909').setTitle('Y-Dimension of detector (in microns).').onChange( update );
	ydimRow.add( new UI.Text( 'Y-Dimension' ).setWidth( '90px' ).setColor( '#666' ) );
	ydimRow.add( ydimIn );
	container.add( ydimRow );
	
	// Detector Thickness
	var detectorRow = new UI.Panel();
	var detectorIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('150').setTitle('Thickness of detector (in microns)').onChange( update );
	detectorRow.add( new UI.Text( 'Detector Thickness' ).setWidth( '90px' ).setColor( '#666' ) );
	detectorRow.add( detectorIn );
	container.add( detectorRow );
	
	// Column Radius
	var colRadRow = new UI.Panel();
	var colRadIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('5.1').setTitle('Column Radius (in microns)').onChange( update );
	colRadRow.add( new UI.Text( 'Column Radius' ).setWidth( '90px' ).setColor( '#666' ) );
	colRadRow.add( colRadIn );
	container.add( colRadRow );
	
	// Column Refractive Index
	var criRow = new UI.Panel();
	var criIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('1.8').setTitle('Refractive index of column material').onChange( update );
	criRow.add( new UI.Text( 'Column Refractive Index' ).setWidth( '90px' ).setColor( '#666' ) );
	criRow.add( criIn );
	container.add( criRow );
	
	// Inter-Columnar Refractive Index
	var icriRow = new UI.Panel();
	var icriIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('1').setTitle('Refractive index of inter-columnar material').onChange( update );
	icriRow.add( new UI.Text( 'Inter-Columnar Refractive Index' ).setWidth( '90px' ).setColor( '#666' ) );
	icriRow.add( icriIn );
	container.add( icriRow );
	
	// Top Surface Absorption Fraction
	var tsafRow = new UI.Panel();
	var tsafIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0.1').setTitle('This parameter dictates how many optical photons are absorbed at the top surface. A fraction oif 1 means every photon that hits the top surface will be absorbed.').onChange( update );
	tsafRow.add( new UI.Text( 'Top Surface Absorption Fraction' ).setWidth( '90px' ).setColor( '#666' ) );
	tsafRow.add( tsafIn );
	container.add( tsafRow );
	
	// Bulk Absorption Coefficient
	var bacRow = new UI.Panel();
	var bacIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0.0001').setTitle('Bulk absorption coefficient (in 1/microns)').onChange( update );
	bacRow.add( new UI.Text( 'Bulk Absorption Coefficient' ).setWidth( '90px' ).setColor( '#666' ) );
	bacRow.add( bacIn );
	container.add( bacRow );
	
	// Surface Roughness Coefficient
	var srcRow = new UI.Panel();
	var srcIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0.2').setTitle('This parameter dictates the roughness of the column walls. A coefficient of 0 means perfectly smooth walls.').onChange( update );
	srcRow.add( new UI.Text( 'Surface Roughness Coefficient' ).setWidth( '90px' ).setColor( '#666' ) );
	srcRow.add( srcIn );
	container.add( srcRow );
	
	// Minimum Distance Next Column
	var mindncRow = new UI.Panel();
	var mindncIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('1').setTitle('Using Columnar Crosstalk (CCT), this is the minimum distance at which next column may be detected (in Microns). The distance to next column is sampled uniformly between the minimum and maximum distance parameters.').onChange( update );
	mindncRow.add( new UI.Text( 'Minimum Distance Next Column' ).setWidth( '90px' ).setColor( '#666' ) );
	mindncRow.add( mindncIn );
	container.add( mindncRow );
	
	// Maximum Distance Next Column
	var maxdncRow = new UI.Panel();
	var maxdncIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('280').setTitle('Using Columnar Crosstalk (CCT), this is the maximum distance at which next column may be detected (in Microns). The distance to next column is sampled uniformly between the minimum and maximum distance parameters.').onChange( update );
	maxdncRow.add( new UI.Text( 'Maximum Distance Next Column' ).setWidth( '90px' ).setColor( '#666' ) );
	maxdncRow.add( maxdncIn );
	container.add( maxdncRow );
	
	// PRF Image X Lower Bound
	var prfxlbRow = new UI.Panel();
	var prfxlbIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0').setTitle('In microns, start of lower bound of X').onChange( update );
	prfxlbRow.add( new UI.Text( 'PRF Image X Lower Bound' ).setWidth( '90px' ).setColor( '#666' ) );
	prfxlbRow.add( prfxlbIn );
	container.add( prfxlbRow );
	
	// PRF Image Y Lower Bound
	var prfylbRow = new UI.Panel();
	var prfylbIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0').setTitle('Y lower bound of PRF image').onChange( update );
	prfylbRow.add( new UI.Text( 'PRF Image Y Lower Bound' ).setWidth( '90px' ).setColor( '#666' ) );
	prfylbRow.add( prfylbIn );
	container.add( prfylbRow );
	
	// PRF Image X Upper Bound
	var prfxubRow = new UI.Panel();
	var prfxubIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('909').setTitle('X upper bound of PRF image').onChange( update );
	prfxubRow.add( new UI.Text( 'PRF Image X Upper Bound' ).setWidth( '90px' ).setColor( '#666' ) );
	prfxubRow.add( prfxubIn );
	container.add( prfxubRow );
	
	// PRF Image Y Upper Bound
	var prfyubRow = new UI.Panel();
	var prfyubIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('909').setTitle('Y upper bound of PRF image').onChange( update );
	prfyubRow.add( new UI.Text( 'PRF Image Y Upper Bound' ).setWidth( '90px' ).setColor( '#666' ) );
	prfyubRow.add( prfyubIn );
	container.add( prfyubRow );
	
	// Light Yield
	var lightRow = new UI.Panel();
	var lightIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0.055').setTitle('light yield (/eV)').onChange( update );
	lightRow.add( new UI.Text( 'Light Yield' ).setWidth( '90px' ).setColor( '#666' ) );
	lightRow.add( lightIn );
	container.add( lightRow );
	
	// Pixel Pitch
	var pixelRow = new UI.Panel();
	var pixelIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('9').setTitle('pixel pitch (in microns)  (max. pixels allowed in PRF image are 501x501. calculate this by {upper bound - lower bound}/pixel pitch.)').onChange( update );
	pixelRow.add( new UI.Text( 'Pixel Pitch' ).setWidth( '90px' ).setColor( '#666' ) );
	pixelRow.add( pixelIn );
	container.add( pixelRow );
	
	// Non-Ideal Sensor Reflectivity
	var nisrRow = new UI.Panel();
	var nisrIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('0.25').setTitle('Fraction of photons to be reflected back').onChange( update );
	nisrRow.add( new UI.Text( 'Non-Ideal Sensor Reflectivity' ).setWidth( '90px' ).setColor( '#666' ) );
	nisrRow.add( nisrIn );
	container.add( nisrRow );
	
	// Machine Number
	var machineRow = new UI.Panel();
	var machineIn = new UI.Input().setWidth( '150px' ).setColor( '#444' ).setFontSize( '12px' ).setValue('1').setTitle('machine number').onChange( update );
	machineRow.add( new UI.Text( 'Machine Number' ).setWidth( '90px' ).setColor( '#666' ) );
	machineRow.add( machineIn );
	container.add( machineRow );
	
	//async code that gets the GPU information when WebMANTIS first loads
	outliner.setValue("Getting GPU information on current node....");
	$.ajax(	//return number of GPUs on node
	{
		type: "POST",
		url: "getGPUInfo.php",
		async: true,
		success: function(data)
		{
			var res;
			var count = 0;
			var response = jQuery.parseJSON(data);
			if(response.ret == "True")
			{
				res = response.text;
				outliner.setValue(res);
			}
		}
	});
	
	//checks the cookie to set the current user correctly
	checkCookie();
	
	$.ajax({ //gets all zipped files currently and appends them to download selector
		type: "POST",
		url: "getZippedFiles.php",
		data: { folder : currentUser},
		async: false,
		success: function(data)
		{
			var response = jQuery.parseJSON(data);
			$.each(response, function(key, val) 
			{
				zipFiles.addOption(val);
			});
		}
	});
	
	$.ajax(
	{ //cleans up any data - Han
		type: "POST",
		url: "cleanUpScript.php",
		data : {folder : currentUser},
	});
			
	function updateUser()
	{
		docCookies.setItem("username", users.getValue(), 31536e3);
		currentUser = users.getValue();
		location.reload();
	}
	
	//check if user has already used previously and sets the
	//corresponding user information
	function checkCookie()
	{
		//var username=getCookie("username");
		var usern = docCookies.getItem("username");
		if (usern!= null && usern != "")
		{
			users.setValue(usern);
			currentUser = usern;
			
			$.ajax(
			{
				type: "POST",
				url: "getInputParameters.php",
				async: false,
				data: 
				{ 
					folder : currentUser
				},
				success: function(data)
				{
					var response = jQuery.parseJSON(data);
					if(response.ret == "True")
					{
						xrayIn.setValue(response.xray);
						prfxubIn.setValue(response.prfxub);
						prfyubIn.setValue(response.prfyub);
						photonsIn.setValue(response.photons);
						mindncIn.setValue(response.mindnc);
						maxdncIn.setValue(response.maxdnc); 
						mindIn.setValue(response.minD);
						maxdIn.setValue(response.maxD);
						numBinsIn.setValue(response.numBins);
						xdimIn.setValue(response.xDim);
						ydimIn.setValue(response.yDim);
						detectorIn.setValue(response.detector);
						colRadIn.setValue(response.colRad);
						criIn.setValue(response.cri);
						icriIn.setValue(response.icri);
						tsafIn.setValue(response.tsaf);
						bacIn.setValue(response.bac);
						prfxlbIn.setValue(response.prfxlb);
						prfylbIn.setValue(response.prfylb);
						lightIn.setValue(response.light);
						pixelIn.setValue(response.pixel);
						nisrIn.setValue(response.nisr);
						srcIn.setValue(response.src);
					}
				}
			});
		}
		else
		{
			users.setValue("Guest");
			currentUser = "Guest";
		}
	}

	function updateZip()
	{
		//TODO
	}
	
	// function that updates the user inputs
	function update() 
	{
		xray2 = xrayIn.getValue();
		prfxub2 = prfxubIn.getValue();
		prfyub2 = prfyubIn.getValue();
		photons2 = photonsIn.getValue();
		mindnc2 = mindncIn.getValue();
		maxdnc2 = maxdncIn.getValue(); 
		mind2 = mindIn.getValue();
		maxd2 = maxdIn.getValue();
		numBins2 = numBinsIn.getValue();
		xdim2 = xdimIn.getValue();
		ydim2 = ydimIn.getValue();
		detector2 = detectorIn.getValue();
		colRad2 = colRadIn.getValue();
		cri2 = criIn.getValue();
		icri2 = icriIn.getValue();
		tsaf2 = tsafIn.getValue();
		bac2 = bacIn.getValue();
		prfxlb2 = prfxlbIn.getValue();
		prfylb2 = prfylbIn.getValue();
		light2 = lightIn.getValue();
		pixel2 = pixelIn.getValue();
		nisr2 = nisrIn.getValue();
		src2 = srcIn.getValue();
	}
	
	//function that updates the progress bar
	var updateProgress = function ()
	{
		$.ajax(
		{
			type: "POST",
			url: "progressBar.php",
			data: { folder : currentUser},
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				if(response.ret == "True")
				{
					percentNow.setValue('Current History: '+response.num);
					currentXray = response.num;
				}
			}
		});
		
		if(currentXray < xray2)
		{
			setTimeout(updateProgress, 500);
		}
		else
		{
			//alert("updateProgress finished");
		}
	}

	// updates the output textbox with latest results
	var updateTextbox = function ()
	{
		$.ajax(
		{
			type: "POST",
			url: "getOutput.php",
			data: { folder : currentUser},
			success: function(data)
			{
				var res;
				var response = jQuery.parseJSON(data);
				
				if(response.ret == "True")
				{
					res = response.text;
					outliner.setValue(res);
				}
			}
		});
		
		if(currentXray < xray2)
		{
			setTimeout(updateTextbox, 500);
		}
		else
		{
			//alert("updateTextbox finished");
		}
	}

	//fetches the data from the backend server
	var fetchData = function ()
	{
		$.ajax(
		{
			type: "POST",
			url: "parseData.php",
			data: { n : fetchCounter, folder : currentUser },
			async: false,
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				
				if(response.ret == "True")
				{
					$.ajax(
					{
						dataType: "json",
						url: response.fileN,
						data: data,
						async: false,
						success: function(data)
						{
							globalJSONData = data;
							editor.addCylinder();
							editor.addSpheres();
							//alert(JSON.stringify(globalJSONData, undefined, 1));
						}
					});
					
					fetchCounter += 1;
				}
			}
		});
		
		if(fetchCounter < photons2)
		{
			setTimeout(fetchData, 500);
		}
		else
		{
			//alert("fetchData finished");
			fetchCounter -= 1;
		}
	}
	
	// reads and shows the PHS image
	var processDetectPNG = function()
	{
		$.ajax(
		{
			type: "POST",
			url: "processDetectPNG.php",
			data: { n : phsCounter, folder : currentUser },
			async: false,
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				if(response.ret == "True")
				{
					var phsimage = currentUser+"/images/PHS_detect"+phsCounter+".png";
					//exist
					var canvas = document.getElementById("phs");
					var context = canvas.getContext("2d");
					var cat = new Image();
					cat.src = phsimage;
					canvas.width = 368;
					canvas.height = 300;
					cat.onload = function() 
					{
						context.drawImage(cat, 0, 0, canvas.width, canvas.height);
					};
					
					phsCounter++;
				}
			}
		});
		
		if(currentXray < xray2)
		{
			setTimeout(processDetectPNG, 500);
		}
		else
		{
			phsCounter = 1;
			var phsimage = currentUser+"/images/PHS_detect"+phsCounter+".png";
			$.get(phsimage)
			.done(function() 
			{	
				//exist
				var canvas = document.getElementById("phs");
				var context = canvas.getContext("2d");
				var cat = new Image();
				cat.src = phsimage;
				canvas.width = 368;
				canvas.height = 300;
				cat.onload = function() 
				{
					context.drawImage(cat, 0, 0, canvas.width, canvas.height);
				};
			});
				
			//alert("processDetectPNG finished");
		}
	}

	// reads and shows the PRF image
	var processMyImagePNG = function ()
	{
		$.ajax(
		{
			type: "POST",
			url: "processMyImagePNG.php",
			data: { n : prfCounter, folder : currentUser },
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				if(response.ret == "True")
				{
					var prfimage = currentUser+"/images/PRF_image"+prfCounter+".png";
					
					//exist
					var canvas = document.getElementById("prf");
					var context = canvas.getContext("2d");
					var cat = new Image();
					cat.src = prfimage;
					canvas.width = 368;
					canvas.height = 300;
					cat.onload = function() 
					{
						context.drawImage(cat, 0, 0, canvas.width, canvas.height);
					};

					
					
					prfCounter++;
				}
			}
		});
		
		if(currentXray < xray2)
		{
			setTimeout(processMyImagePNG, 500);
		}
		else
		{
			prfCounter = 1;
			var prfimage = currentUser+"/images/PRF_image"+prfCounter+".png";
			$.get(prfimage)
			.done(function() 
			{	
				//exist
				var canvas = document.getElementById("prf");
				var context = canvas.getContext("2d");
				var cat = new Image();
				cat.src = prfimage;
				canvas.width = 368;
				canvas.height = 300;
				cat.onload = function() 
				{
					context.drawImage(cat, 0, 0, canvas.width, canvas.height);
				};
			});
			//alert("processMyImagePNG finished");
		}
	}
	
	// left click on the PHS button
	function phsLeftClick()
	{
		var phsimage = currentUser+"/images/PHS_detect"+(phsCounter-1)+".png";
		$.get(phsimage)
		.done(function() 
		{	
			//exist
			var canvas = document.getElementById("phs");
			var context = canvas.getContext("2d");
			var cat = new Image();
			cat.src = phsimage;
			canvas.width = 368;
			canvas.height = 300;
			cat.onload = function() 
			{
				context.drawImage(cat, 0, 0, canvas.width, canvas.height);
			};
			phsCounter -= 1;
		})
		.fail(function() { 
			//alert("button.lPHS failed");
		});	
	}
	
	// right click on the PHS button
	function phsRightClick()
	{
		var phsimage = currentUser+"/images/PHS_detect"+(phsCounter+1)+".png";
		$.get(phsimage)
		.done(function() 
		{	
			//exist
			var canvas = document.getElementById("phs");
			var context = canvas.getContext("2d");
			var cat = new Image();
			cat.src = phsimage;
			canvas.width = 368;
			canvas.height = 300;
			cat.onload = function() 
			{
				context.drawImage(cat, 0, 0, canvas.width, canvas.height);
			};
			phsCounter += 1;
		})
		.fail(function() { 
			//alert("button.lPHS failed");
		});	
	}
	
	// left click on the PRF button
	function prfLeftClick()
	{
		var prfimage = currentUser+"/images/PRF_image"+(prfCounter-1)+".png";
		$.get(prfimage)
		.done(function() 
		{	
			//exist
			var canvas = document.getElementById("prf");
			var context = canvas.getContext("2d");
			var cat = new Image();
			cat.src = prfimage;
			canvas.width = 368;
			canvas.height = 300;
			cat.onload = function() 
			{
				context.drawImage(cat, 0, 0, canvas.width, canvas.height);
			};
			prfCounter -= 1;
		})
		.fail(function() { 
			//alert("button.lPHS failed");
		});	
	}
	
	// right click on the PRF button
	function prfRightClick()
	{
		var prfimage = currentUser+"/images/PRF_image"+(prfCounter+1)+".png";
		$.get(prfimage)
		.done(function() 
		{	
			//exist
			var canvas = document.getElementById("prf");
			var context = canvas.getContext("2d");
			var cat = new Image();
			cat.src = prfimage;
			canvas.width = 368;
			canvas.height = 300;
			cat.onload = function() 
			{
				context.drawImage(cat, 0, 0, canvas.width, canvas.height);
			};
			prfCounter += 1;
		})
		.fail(function() { 
			//alert("button.lPHS failed");
		});	
	}
	
	// events
	signals.saveInputArguments.add( function ()
	{
		update();
	});
	
	// starts the execution of program
	signals.executeStarted.add( function () 
	{
		mypid.setValue('Process ID: '+myPID);
		updateProgress();
		updateTextbox();
		fetchData();
		processDetectPNG();
		processMyImagePNG();
	});
	
	return container;
}
