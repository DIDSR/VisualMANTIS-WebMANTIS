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
//	Filename:	my_image_cb.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file handles displaying the Point Response Function image.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************
 * myimage_cb
 * 
 * Callback that handles displaying
 * the myimage data file on the
 * GUI.
 * ***************************/
void myimage_cb(void *)
{
	FILE * imagefp, *console;
	char filename[100];
	char imagename[100];
	
	sprintf(filename, "myimage%d.dat", imageCBcount);
	sprintf(imagename, "PRF_han%d.png", imageCBcount);

	//if myimage.dat exists and simulation has started
	imagefp = fopen(filename, "r");
	if(imagefp != NULL && pressed == 1)
	{
		//image data written
		fclose(imagefp);
		
		//uses popen to open a pipe to console for gnuplot to generate png image
		console = popen("gnuplot", "w");
		fprintf(console, "set terminal png enh size 500, 400 crop\n");
		fprintf(console, "set out \'%s\'\n", imagename);
		fprintf(console, "set pm3d map\n");
		fprintf(console, "set log cb\n");
		fprintf(console, "set cbrange[1e-3:1e2]\n");
		fprintf(console, "set log z\n");
		fprintf(console, "set si sq\n");
		fprintf(console, "set cbtics format \'1e%%T\'\n");
		fprintf(console, "set xlabel \'x label\'\n");
		fprintf(console, "set ylabel \'y label\'\n");
		fprintf(console, "unset key\n");
		fprintf(console, "set palette rgbformulae 33,13,10\n");
		fprintf(console, "splot \'%s\'\n", filename);
		fprintf(console, "exit\n");
		fclose(console);
		//printf("finished writing: %s\n", imagename);
		
		//gets png image to display
		image = new Fl_PNG_Image(imagename);
		image = (Fl_PNG_Image*) image->copy(b3->w()-100, b3->h()-200);
		
		b3->position(1150, 600);
		b3->box(FL_NO_BOX);
		b3->image(image);
		b3->redraw();

		//updates label value
		imageoutput->value(imagename);
		imageoutput->redraw();
		imageCBcount ++;
	}
	left->redraw();
	right->redraw();
	
	//recursive call every 1 second
	Fl::repeat_timeout(1.0, myimage_cb);
}

/******************************
 * left_cb
 * 
 * Callback that handles
 * the left button click to
 * change the PRF image.
 * ***************************/
void left_cb(void *)
{
	FILE * imagefp;

	imageLRCBcount --;
	char imagename[100];
	sprintf(imagename, "PRF_han%d.png", imageLRCBcount);

	imagefp = fopen(imagename, "r");
	if(imagefp != NULL)
	{
		image = new Fl_PNG_Image(imagename);
		image = (Fl_PNG_Image*) image->copy(500, 350);
		b3->image(image);
		b3->redraw();
		imageoutput->value(imagename);
		//imageoutput->redraw();
		fclose(imagefp);
	}
	else
	{
		imageLRCBcount ++;
	}
	left->redraw();
	right->redraw();
}

/******************************
 * right_cb
 * 
 * Callback that handles
 * the right button click to
 * change the PRF image.
 * ***************************/
void right_cb(void *)
{
	FILE * imagefp;

	imageLRCBcount ++;
	char imagename[100];
	sprintf(imagename, "PRF_han%d.png", imageLRCBcount);

	imagefp = fopen(imagename, "r");
	if(imagefp != NULL && pressed == 1)
	{
		image = new Fl_PNG_Image(imagename);
		image = (Fl_PNG_Image*) image->copy(500, 350);
		b3->image(image);
		b3->redraw();
		imageoutput->value(imagename);
		//imageoutput->redraw();
		fclose(imagefp);
	}
	else
	{
		imageLRCBcount --;
	}
	left->redraw();
	right->redraw();
}
