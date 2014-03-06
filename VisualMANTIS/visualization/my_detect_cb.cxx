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
//	Filename:	my_detect_cb.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file handles displaying the Pulse Height Spectrum images.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/******************************
 * detect_cb
 * 
 * Callback that handles displaying
 * the detect data file on the
 * GUI.
 * ***************************/
void detect_cb(void *)
{
	FILE * detectfp, *console;
	char filename[100];
	char detectname[100];
	
	sprintf(filename, "detect%d.dat", detectCBcount);
	sprintf(detectname, "PHS_detect%d.png", detectCBcount);
	
	//if detect.dat exists and simulation has started
	detectfp = fopen(filename, "r");
	if(detectfp != NULL && pressed == 1)
	{
		//image data written
		fclose(detectfp);
		
		//uses popen to open a pipe to console for gnuplot to generate png image
		console = popen("gnuplot", "w");
		fprintf(console, "set terminal png enh size 500, 500 crop\n");
		fprintf(console, "set out \'%s\'\n", detectname);
		fprintf(console, "B=%d\n", (atoi(maxD->value()) - atoi(minD->value())) / atoi(bins->value()));
		fprintf(console, "set style line 1 linetype 1 linewidth 4\n");
		fprintf(console, "set xlabel \'# detected photons per primary\'\n");
		fprintf(console, "set ylabel \'Frequency\'\n");
		fprintf(console, "set yrange[0:30000]\n");
		fprintf(console, "set si sq\n");
		fprintf(console, "unset key\n");
		fprintf(console, "plot \'%s\' u ($0*B):1 w l ls 1\n", filename);
		fprintf(console, "exit\n");
		fclose(console);
		//printf("finished writing: %s\n", detectname);
		
		//gets png image to display
		image2 = new Fl_PNG_Image(detectname);
		image2 = (Fl_PNG_Image*) image2->copy(b2->w()-100, b2->h()-200);
		//printf("b2->w: %d b2->h: %d\n", b2->w(), b2->h());
		
		b2->position(1150, 100);
		b2->box(FL_NO_BOX);
		b2->image(image2);
		b2->redraw();
		
		//updates label value		
		detectOutput->value(detectname);
		detectOutput->redraw();
		detectCBcount ++;
	}
	Fl::repeat_timeout(1.0, detect_cb);
}

/******************************
 * detectRight_cb
 * 
 * Callback that handles
 * the right button click to
 * change the detect image.
 * ***************************/
void detectRight_cb(void *)
{
	FILE * detectfp;;

	detectLRCBcount ++;
	char imagename[100];
	sprintf(imagename, "PHS_detect%d.png", detectLRCBcount);

	//if the image exists
	detectfp = fopen(imagename, "r");
	if(detectfp != NULL)
	{
		image2 = new Fl_PNG_Image(imagename);
		image2 = (Fl_PNG_Image*) image2->copy(500, 350);
		b2->image(image2);
		b2->redraw();
		detectOutput->value(imagename);
		//imageoutput->redraw();
		fclose(detectfp);
	}
	else
	{
		detectLRCBcount --;
	}
	detectLeft->redraw();
	detectRight->redraw();
}

/******************************
 * detectLeft_cb
 * 
 * Callback that handles
 * the left button click to
 * change the detect image.
 * ***************************/
void detectLeft_cb(void *)
{
	FILE * detectfp;;

	detectLRCBcount --;
	char imagename[100];
	sprintf(imagename, "PHS_detect%d.png", detectLRCBcount);

	//if the image exists
	detectfp = fopen(imagename, "r");
	if(detectfp != NULL)
	{
		image2 = new Fl_PNG_Image(imagename);
		image2 = (Fl_PNG_Image*) image2->copy(500, 350);
		b2->image(image2);
		b2->redraw();
		detectOutput->value(imagename);
		//imageoutput->redraw();
		fclose(detectfp);
	}
	else
	{
		detectLRCBcount ++;
	}
	detectLeft->redraw();
	detectRight->redraw();
}
