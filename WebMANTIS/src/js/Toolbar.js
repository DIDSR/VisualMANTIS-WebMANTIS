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
//	Filename:	Toolbar.php
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Toolbar is at the bottom of the screen and controls the visuals in the main screen
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
var Toolbar = function ( editor ) {

	var signals = editor.signals;

	var container = new UI.Panel();
	container.setPosition( 'absolute' );
	container.setClass( 'toolbar' );

	var buttons = new UI.Panel();
	container.add( buttons );

	//left and right button clicks changes the histories
	var left = new UI.Button( '<--' ).onClick( function () 
	{
		leftClick();
	});
	buttons.add( left );
	var visualCounter = new UI.Text( 'Visualization: ' )
	buttons.add( visualCounter );
	var right = new UI.Button( '-->' ).onClick( function () 
	{
		rightClick();
	});
	buttons.add( right );
	
	// alters the opacity of the cylinders
	var opacity = new UI.Number().setRange( 0, 1 ).onChange( updateOpacity );
	opacity.dom.style.width = '42px';
	buttons.add( new UI.Text( 'Opacity: ' ) );
	buttons.add( opacity );
	
	// photon history - slider that can individually show the photons
	// trajectory in single timesteps
	var historyCounter = new UI.Text( 'Photon History: 0 ' )
	var pHistory = new UI.Slider(historyCounter).setRange(0, 1).setStep(1).setValue(0).onChange(updateHistory);
	buttons.add(historyCounter);
	buttons.add(new UI.Text());
	buttons.add( pHistory );
	
	//function for left click above
	function leftClick()
	{
		if((fetchCounter-1) >= 0 && (fetchCounter-1) < photons2)
		{
			//hide prior stuff
			for(var j=0;j<objects[fetchCounter].length;j++)
			{
				objects[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<spheres[fetchCounter].length;j++)
			{
				spheres[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<lines[fetchCounter].length;j++)
			{
				lines[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<ceilings[fetchCounter].length;j++)
			{
				ceilings[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<floors[fetchCounter].length;j++)
			{
				floors[fetchCounter][j].visible = false;
			}
		
			fetchCounter -= 1;
			//show new stuff
			for(var j=0;j<objects[fetchCounter].length;j++)
			{
				objects[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<spheres[fetchCounter].length;j++)
			{
				spheres[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<lines[fetchCounter].length;j++)
			{
				lines[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<ceilings[fetchCounter].length;j++)
			{
				ceilings[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<floors[fetchCounter].length;j++)
			{
				floors[fetchCounter][j].visible = true;
			}
			signals.nextHistory.dispatch();
			signals.renderScene.dispatch();
		}
	}
	
	//function for right click above
	function rightClick()
	{
		if((fetchCounter+1) >= 0 && (fetchCounter+1) < photons2)
		{
			//hide prior stuff
			for(var j=0;j<objects[fetchCounter].length;j++)
			{
				objects[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<spheres[fetchCounter].length;j++)
			{
				spheres[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<lines[fetchCounter].length;j++)
			{
				lines[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<ceilings[fetchCounter].length;j++)
			{
				ceilings[fetchCounter][j].visible = false;
			}
			
			for(var j=0;j<floors[fetchCounter].length;j++)
			{
				floors[fetchCounter][j].visible = false;
			}
		
			fetchCounter += 1;
			//show new stuff
			for(var j=0;j<objects[fetchCounter].length;j++)
			{
				objects[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<spheres[fetchCounter].length;j++)
			{
				spheres[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<lines[fetchCounter].length;j++)
			{
				lines[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<ceilings[fetchCounter].length;j++)
			{
				ceilings[fetchCounter][j].visible = true;
			}
			
			for(var j=0;j<floors[fetchCounter].length;j++)
			{
				floors[fetchCounter][j].visible = true;
			}
			signals.nextHistory.dispatch();
			signals.renderScene.dispatch();
		}
	}
	
	//changes opacity
	function updateOpacity() 
	{
		for(var j=0;j<objects[fetchCounter].length;j++)
		{
			objects[fetchCounter][j].material.opacity = opacity.getValue();
		}
		
		signals.renderScene.dispatch();
	}

	// updates the path of each individual photon
	function updateHistory() 
	{
		//sets all previous cylinders to be visible up to the current one
		for(var j=0;j<objects[fetchCounter].length;j++)
		{
			if(j <= pHistory.getValue())
			{
				objects[fetchCounter][j].visible = true;
			}
			else
			{
				objects[fetchCounter][j].visible = false;
			}
		}
		
		//sets all previous photons to be visible up to the current one
		for(var j=0;j<spheres[fetchCounter].length;j++)
		{
			if(j <= pHistory.getValue())
			{
				spheres[fetchCounter][j].visible = true;
			}
			else
			{
				spheres[fetchCounter][j].visible = false;
			}
		}
		
		//sets all previous lines to be visible up to the current one
		for(var j=0;j<lines[fetchCounter].length;j++)
		{
			if(j <= pHistory.getValue())
			{
				lines[fetchCounter][j].visible = true;
			}
			else
			{
				lines[fetchCounter][j].visible = false;
			}
		}
		
		//sets all previous floors and ceilings to be visible up to the
		//current one
		for(var j=0;j<ceilings[fetchCounter].length;j++)
		{
			if(j <= pHistory.getValue())
			{
				ceilings[fetchCounter][j].visible = true;
			}
			else
			{
				ceilings[fetchCounter][j].visible = false;
			}
		}
		
		for(var j=0;j<floors[fetchCounter].length;j++)
		{
			if(j <= pHistory.getValue())
			{
				floors[fetchCounter][j].visible = true;
			}
			else
			{
				floors[fetchCounter][j].visible = false;
			}
		}
		
		signals.renderScene.dispatch();
	}
	
	signals.nextHistory.add( function () 
	{
		visualCounter.setValue("Visualization: "+fetchCounter+"  ");
		pHistory.setRange(0, spheres[fetchCounter].length);
	});

	return container;

}
