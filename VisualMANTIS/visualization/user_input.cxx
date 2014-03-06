///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        visualMANTIS v1.0		     //
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
//	Filename:	user_input.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file handles geetting inputs from the user.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************
 * save_cb()
 * 
 * Callback corresponding to
 * "Save" button that saves
 * user input changes.
 * ***************************/
void save_cb(Fl_Widget *)
{
	int start = 1;
	char alertline[1000];
	
	//checks if the following variables exceed the limit and alerts the user.
	//sets start to 0 if there are errors
	/*if(atof(maxD->value()) > atof(xRay->value()) * atof(lightYield->value()))
	{
		sprintf(alertline, "Maximum number of photons that can be detected must be less than %lf (N * lightYield)\n", atof(xRay->value()) * atof(lightYield->value()));
		start = 0;
	}*/
	if(atoi(bins->value()) > 1000)
	{
		strcat(alertline, "Number of bins must be less than 1000\n");
		start = 0;
	}
	if(atoi(pixelPitch->value()) > 501)
	{
		strcat(alertline, "Pixel pitch must be less than 501\n");
		start = 0;
	}
	
	//if no issues above, saves modified inputs
	if(start == 1)
	{
		FILE * tempFP;
	
		tempFP = fopen("temp.txt", "w");
	
		fprintf(tempFP, "%s\n", xRay->value());
		fprintf(tempFP, "%s\n", minD->value());
		fprintf(tempFP, "%s\n", maxD->value());
		fprintf(tempFP, "%s\n", bins->value());
		fprintf(tempFP, "%s\n", xDim->value());
		fprintf(tempFP, "%s\n", yDim->value());
		fprintf(tempFP, "%s\n", dThickness->value());
		fprintf(tempFP, "%s\n", colRad->value());
		fprintf(tempFP, "%s\n", colRefrac->value());
		fprintf(tempFP, "%s\n", interColRefrac->value());
		fprintf(tempFP, "%s\n", absorpFrac->value());
		fprintf(tempFP, "%s\n", bulkAbsorpCoeff->value());
		fprintf(tempFP, "%s\n", surRoughCoeff->value());
		fprintf(tempFP, "%s\n", minDistance->value());
		fprintf(tempFP, "%s\n", maxDistance->value());
		fprintf(tempFP, "%s\n", xLower->value());
		fprintf(tempFP, "%s\n", yLower->value());
		fprintf(tempFP, "%s\n", xUpper->value());
		fprintf(tempFP, "%s\n", yUpper->value());
		fprintf(tempFP, "%s\n", lightYield->value());
		fprintf(tempFP, "%s\n", pixelPitch->value());
		fprintf(tempFP, "%s\n", niSensor->value());
		fprintf(tempFP, "%s\n", flag->value());
		fprintf(tempFP, "%s\n", machine->value());
		fprintf(tempFP, "%s\n", numPhotonHistories->value());
		
		fclose(tempFP);
		
		subWindow->hide();
	}
	else
	{
		fl_alert(alertline, 10.0);
	}
}

/******************************
 * cancel_cb()
 * 
 * Callback for "Cancel" button
 * in input subWindow. Does not 
 * save modifications in user 
 * inputs.
 * ***************************/
void cancel_cb(Fl_Widget *)
{
	FILE * tempFP;
	char line[100];
	
	//reads the original inputs from temp.txt
	tempFP = fopen("temp.txt", "r");
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	xRay->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	minD->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	maxD->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	bins->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	xDim->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	yDim->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	dThickness->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	colRad->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	colRefrac->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	interColRefrac->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	absorpFrac->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	bulkAbsorpCoeff->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	surRoughCoeff->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	minDistance->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	maxDistance->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	xLower->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	yLower->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	xUpper->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	yUpper->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	lightYield->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	pixelPitch->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	niSensor->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	flag->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	machine->value(line);
	
	fgets(line, 100, tempFP);
	line[strlen(line)-1]='\0';	
	numPhotonHistories->value(line);
	
	fclose(tempFP);
	
	subWindow->hide();
}

/******************************
 * hide_cb()
 * 
 * Callback when user clicks
 * button to show additional
 * inputs in input subwindow
 * ***************************/
void hide_cb(Fl_Widget *) 
{
	//shows invisible inputs
	if(hideCounter == 0)
	{
		subWindow->resize(0, 0, 475, 810);
		saveButton->resize(5, 780, 80, 25);
		cancelButton->resize(90, 780, 80, 25);
		colRefrac->show();	
		interColRefrac->show();	
		absorpFrac->show();	
		bulkAbsorpCoeff->show();	
		surRoughCoeff->show();	
		minDistance->show();	
		maxDistance->show();	
		xLower->show();
		yLower->show();
		xUpper->show();
		yUpper->show();
		lightYield->show();
		pixelPitch->show();	
		niSensor->show();
		flag->show();
		machine->show();
		numPhotonHistories->show();
		hideCounter = 1;
	}
	//hides inputs
	else if(hideCounter == 1)
	{
		subWindow->resize(0, 0, 475, 300);
		saveButton->resize(5, 270, 80, 25);
		cancelButton->resize(90, 270, 80, 25);
		colRefrac->hide();	
		interColRefrac->hide();	
		absorpFrac->hide();	
		bulkAbsorpCoeff->hide();	
		surRoughCoeff->hide();	
		minDistance->hide();	
		maxDistance->hide();	
		xLower->hide();
		yLower->hide();
		xUpper->hide();
		yUpper->hide();
		lightYield->hide();
		pixelPitch->hide();	
		niSensor->hide();
		flag->hide();
		machine->hide();
		numPhotonHistories->hide();
		hideCounter = 0;
	}
}

/**************************************************************
 * initParameterWindow()
 * 
 * Initializes and creates the subwindow for user inputs.
 * ***********************************************************/
void initParameterWindow()
{
	//Window for input parameters
	subWindow = new Fl_Window(0, 0, 475, 300, "Parameters");
	
	//visible inputs 
		xRay = new Fl_Input(295, 30, 165, 25, "X-Ray Histories:");
      	xRay->tooltip("Number of x-ray histories to be simulated (N)");
      	xRay->static_value("100000");
  
		minD = new Fl_Input(295, 60, 165, 25, "Min Detect:");
      	minD->tooltip("Minimum number of optical photons that can be detected");
      	minD->static_value("0");
 
		maxD = new Fl_Input(295, 90, 165, 25, "Max Detect:");
     	maxD->tooltip("Maximum (N*yield) number of optical photons that can be detected");
      	maxD->static_value("1400");
 
		bins = new Fl_Input(295, 120, 165, 25, "Number of Bins:");
      	bins->tooltip("Number of bins for storing pulse height specturm (maximum value = 1000)");
      	bins->static_value("140");
	 
		xDim = new Fl_Input(295, 150, 165, 25, "X-Dimension:");
      	xDim->tooltip("X-dimension of detector (in microns)");
      	xDim->static_value("909.0");
  
		yDim = new Fl_Input(295, 180, 165, 25, "Y-Dimension:");
      	yDim->tooltip("Y-dimension of detector (in microns)");
      	yDim->static_value("909.0");

		dThickness = new Fl_Input(295, 210, 165, 25, "Detector Thickness:");
      	dThickness->tooltip("Thickness of detector (in microns)");
      	dThickness->static_value("150.0");
 
		colRad = new Fl_Input(295, 240, 165, 25, "Column Radius:");
      	colRad->tooltip("Column radius (in microns)");
      	colRad->static_value("5.1");
		
		hideButton = new Fl_Button(5, 240, 165, 25, "Hide/Show More");
    	hideButton->callback((Fl_Callback*) hide_cb);
    	
    	saveButton = new Fl_Button(5, 270, 80, 25, "Save");
    	saveButton->callback((Fl_Callback*) save_cb);
    	
    	cancelButton = new Fl_Button(90, 270, 80, 25, "Cancel");
    	cancelButton->callback((Fl_Callback*) cancel_cb);
    	
	//hidden inputs
		colRefrac = new Fl_Input(295, 270, 165, 25, "Column Refractive Index:");
      	colRefrac->tooltip("Refractive index of column material");
      	colRefrac->static_value("1.8");
		colRefrac->hide();
		
		interColRefrac = new Fl_Input(295, 300, 165, 25, "Inter-Columnar Refractive Index:");
      	interColRefrac->tooltip("Refractive index of inter-columnar material");
      	interColRefrac->static_value("1.0");
		interColRefrac->hide();
		
		absorpFrac = new Fl_Input(295, 330, 165, 25, "Top Surface Absorption Fraction:");
      	absorpFrac->tooltip("Top surface absorption fraction");
      	absorpFrac->static_value("0.1");
		absorpFrac->hide();
		
		bulkAbsorpCoeff = new Fl_Input(295, 360, 165, 25, "Bulk Absorption Coefficient:");
      	bulkAbsorpCoeff->tooltip("Bulk absorption coefficient (in 1/microns)");
      	bulkAbsorpCoeff->static_value("1e-4");
		bulkAbsorpCoeff->hide();
		
		surRoughCoeff = new Fl_Input(295, 390, 165, 25, "Surface Roughness Coefficient:");
      	surRoughCoeff->tooltip("Surface roughness coefficient");
      	surRoughCoeff->static_value("0.2");
		surRoughCoeff->hide();
		
		minDistance = new Fl_Input(295, 420, 165, 25, "Minimum Distance Next Column:");
      	minDistance->tooltip("Minimum distance to the next column (in microns)");
      	minDistance->static_value("1.0");
		minDistance->hide();
		
		maxDistance = new Fl_Input(295, 450, 165, 25, "Maximum Distance to Next Column:");
      	maxDistance->tooltip("Maximum distance to the next column (in microns)");
	    maxDistance->static_value("280.0");
		maxDistance->hide();
		
		xLower = new Fl_Input(295, 480, 165, 25, "PRF Image X Lower Bound");
      	xLower->tooltip("X lower bound of PRF image");
      	xLower->static_value("0.0");
		xLower->hide();
		
		yLower = new Fl_Input(295, 510, 165, 25, "PRF Image Y Lower Bound");
      	yLower->tooltip("Y lower bound of PRF image");
      	yLower->static_value("0.0");
		yLower->hide();
		
		xUpper = new Fl_Input(295, 540, 165, 25, "PRF Image X Upper Bound");
      	xUpper->tooltip("X upper bound of PRF image");
      	xUpper->static_value("909.0");
		xUpper->hide();
		
		yUpper = new Fl_Input(295, 570, 165, 25, "PRF Image Y Upper Bound");
      	yUpper->tooltip("Y upper bound of PRF image");
      	yUpper->static_value("909.0");
		yUpper->hide();
		
		lightYield = new Fl_Input(295, 600, 165, 25, "Light Yield:");
      	lightYield->tooltip("Light yield (/eV)");
      	lightYield->static_value("0.055");
		lightYield->hide();
		
		pixelPitch = new Fl_Input(295, 630, 165, 25, "Pixel Pitch:");
      	pixelPitch->tooltip("Pixel pitch (in microns)  (max. pixels allowed in PRF image are 501x501. calculate this by upper bound - lower bound/pixel pitch.)");
	  	pixelPitch->static_value("9");
		pixelPitch->hide();
		
		niSensor = new Fl_Input(295, 660, 165, 25, "Non-Ideal Sensor Reflectivity:");
      	niSensor->tooltip("Non-ideal sensor reflectivity");
      	niSensor->static_value("0.25");
		niSensor->hide();
		
		flag = new Fl_Input(295, 690, 165, 25, "GPU (1) or CPU (0) flag:");
      	flag->tooltip("Flag for running in the GPU (1) or only in the CPU (0)");
      	flag->static_value("1");
		flag->hide();
		
		machine = new Fl_Input(295, 720, 165, 25, "Machine Number:");
      	machine->tooltip("Machine number");
      	machine->static_value("1");
		machine->hide();
		
      	numPhotonHistories = new Fl_Input(295, 750, 165, 25, "Number Photon Histories Visualization");
      	numPhotonHistories->tooltip("Number of photon histories to keep for simulation");
      	numPhotonHistories->static_value("10");
      	numPhotonHistories->hide();
      	
      	//writes data to a temporary file
      	FILE * tempFP;
		tempFP = fopen("temp.txt", "w");
		fprintf(tempFP, "%s\n", xRay->value());
		fprintf(tempFP, "%s\n", minD->value());
		fprintf(tempFP, "%s\n", maxD->value());
		fprintf(tempFP, "%s\n", bins->value());
		fprintf(tempFP, "%s\n", xDim->value());
		fprintf(tempFP, "%s\n", yDim->value());
		fprintf(tempFP, "%s\n", dThickness->value());
		fprintf(tempFP, "%s\n", colRad->value());
		fprintf(tempFP, "%s\n", colRefrac->value());
		fprintf(tempFP, "%s\n", interColRefrac->value());
		fprintf(tempFP, "%s\n", absorpFrac->value());
		fprintf(tempFP, "%s\n", bulkAbsorpCoeff->value());
		fprintf(tempFP, "%s\n", surRoughCoeff->value());
		fprintf(tempFP, "%s\n", minDistance->value());
		fprintf(tempFP, "%s\n", maxDistance->value());
		fprintf(tempFP, "%s\n", xLower->value());
		fprintf(tempFP, "%s\n", yLower->value());
		fprintf(tempFP, "%s\n", xUpper->value());
		fprintf(tempFP, "%s\n", yUpper->value());
		fprintf(tempFP, "%s\n", lightYield->value());
		fprintf(tempFP, "%s\n", pixelPitch->value());
		fprintf(tempFP, "%s\n", niSensor->value());
		fprintf(tempFP, "%s\n", flag->value());
		fprintf(tempFP, "%s\n", machine->value());
		fprintf(tempFP, "%s\n", numPhotonHistories->value());
		fclose(tempFP);
		
	subWindow->end();
	subWindow->show();
	subWindow->hide();
}
