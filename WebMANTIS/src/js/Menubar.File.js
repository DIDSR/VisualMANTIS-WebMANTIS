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
//	Filename:	Menubar.File.js
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Code that forms the File option in the menubar
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
Menubar.File = function ( editor ) {

	var container = new UI.Panel();
	container.setClass( 'menu' );
	container.onClick( function () 
	{
		if(openvar == 0)
		{ 
			options.setDisplay( 'block' );
			openvar = 1;
		}
		else
		{
			options.setDisplay( 'none' );
			openvar = 0;
		}
	});

	var title = new UI.Panel();
	title.setTextContent( 'File' ).setColor( '#666' );
	title.setMargin( '0px' );
	title.setPadding( '8px' );
	container.add( title );

	//

	var options = new UI.Panel();
	options.setClass( 'options' );
	options.setDisplay( 'none' );
	container.add( options );

	// new

	var option = new UI.Panel();
	option.setClass( 'option' );
	option.setTextContent( 'New' );
	option.onClick( function () {

		if ( confirm( 'Are you sure?' ) ) {

			if ( localStorage.threejsEditor !== undefined ) {

				delete localStorage.threejsEditor;

			}

			location.href = location.pathname;

		}

	} );
	options.add( option );
	options.add( new UI.HorizontalRule() );
	
	// execute visualmantis
	var option = new UI.Panel();
	option.setClass( 'option' );
	option.setTextContent( 'Execute' );
	option.onClick( function () 
	{
		editor.saveInputArguments();
		
		$.ajax({	//writes parameters to input files for visualmantis
			type: "POST",
			url: "writeInputParameters.php",
			async: false,
			data: 
			{ 
				xRay : xray2,
				minDetect : mind2,
				maxDetect : maxd2,
				numBins : numBins2, 
				xDim : xdim2,
				yDim : ydim2,
				detThick : detector2,
				colRad : colRad2,
				cri : cri2,
				icri : icri2,
				tsaf : tsaf2,
				bac : bac2,
				src : src2,
				mindnc : mindnc2,
				maxdnc : maxdnc2,
				prfxlb : prfxlb2,
				prfylb : prfylb2,
				prfxub : prfxub2,
				prfyub : prfyub2,
				light : light2,
				pixelPitch : pixel2,
				isr : nisr2,
				histories : photons2,
				folder : currentUser
			}
		});
								
		//editor.executeStart();
		$.ajax(
		{
			type: "POST",
			url: "executeVisualmantis.php",
			data: { folder : currentUser},
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				myPID = response.pid;
				editor.executeStart();
				//alert("Successful");
			}
		});
	});
	options.add( option );
	options.add( new UI.HorizontalRule() );
	
	// stop job 
	var stopj = new UI.Panel();
	stopj.setClass( 'option' );
	stopj.setTextContent( 'Stop' );
	stopj.onClick( function () 
	{
		$.ajax(
		{
			type: "POST",
			url: "killJob.php",
			async: false,
			data: { folder : currentUser},
			success: function(data)
			{
				var response = jQuery.parseJSON(data);
				if(response.ret == '1')
				{
					alert("Job killed");
					//startSimulation = 0;
				}
				else
				{
					alert('No GPU job to kill');
				}
			}
		});
	});
	options.add( stopj );
	options.add( new UI.HorizontalRule() );
	return container;

}
