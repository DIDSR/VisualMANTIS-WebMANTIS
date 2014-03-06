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
//	Filename:	visualizationGUI.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file is the main GUI files that controls all other aspects of the visualization.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************
 * FDA
 * 
 * GUI for hybridMANTIS
 * 
 * 
 * 
 * 
 * **********************************************************/
#include "visualizationGUI.h"

/*  Poor man's approximation of PI */
#define PI 3.1415926535898
/*  Macro for sin & cos in degrees */
#define Cos(th) cos(PI/180.0*(th))
#define Sin(th) sin(PI/180.0*(th))
/*  D degrees of rotation */
#define DEF_D 5

#define UP 		1
#define DOWN 	0
#define X_AXIS	0
#define Y_AXIS	1
#define Z_AXIS	2

//structures used to save photon information
struct cylinders
{
	int size;
	float *x;
	float *y;
	float *z;
};

struct photons
{
	float *px;
	float *py;
	float *pz;
};

struct ceil
{
	int num;
	float *x;
	float *y;
	float *z;
};
struct histStruct
{
	float x[1000];			// x-coordinate
	float y[1000];			// y-coordinate
	float z[1000];			// z-coordinate
	float Xc[1000];			// x coordinate of cylinder
	float Yc[1000];			// y coordinate of cylinder
	int terminated[1000];	//counter to indicate type of termination
	int ceiling[1000];		// 1 or 0 to indicate whether a ceiling should be draw
	
	int numTerminates;
	int dynamicCounter;
	int histCounter;	// counter to keep track of number of histories
};

//default simulation variables
char * def_xRay =  (char *)"100000";
char * def_minD = (char *)"0";
char * def_maxD = (char *)"1400";
char * def_bins = (char *)"140";
char * def_xDim = (char *)"909.0";
char * def_yDim = (char *)"909.0";
char * def_dThickness = (char *)"150.0";
char * def_colRad = (char *)"5.1";
char * def_colRefrac = (char *)"1.8";
char * def_interColRefrac = (char *)"1.0";
char * def_absorpFrac = (char *)"0.1";
char * def_bulkAbsorpCoeff = (char *)"1e-4";
char * def_surRoughCoeff = (char *)"0.2";
char * def_minDistance = (char *)"1.0";
char * def_maxDistance = (char *)"280.0";
char * def_xLower = (char *)"0.0";
char * def_yLower = (char *)"0.0";
char * def_xUpper = (char *)"909.0";
char * def_yUpper = (char *)"909.0";
char * def_lightYield = (char *)"0.055";
char * def_pixelPitch = (char *)"9";
char * def_niSensor = (char *)"0.25";
char * def_flag = (char *)"1";
char * def_machine = (char *)"1";
      	
int pressed = 0;

int numHistStruct = 10;

//GLvoid *font_style = GLUT_BITMAP_TIMES_ROMAN_24;

const float ALPHA = 0.65f; //The opacity of each face of the cylinder

/* vector */
typedef float vec3_t[3]; 

/* amount to rotate about axis */
float rotate = 0.0f;
float rotate2 = 0.0f;

/* vector which describes the axis to rotate about */
vec3_t axis = {1.0, 0.0, 0.0};

int bitmapHeight=35;
//int font=(int) GLUT_BITMAP_TIMES_ROMAN_10;

float vertx, verty, vertz;
int mainWindow;
int mainWindow2;
int mainWindow3;
struct histStruct * photonHistory;
struct histStruct * photonHistory3;

int oldInt = -99;
int newInt = 0;
int renderCount = 1;
int glDataCounter = 0;
int glDataCounter2 = 0;
int glDataCounter3 = 0;

//3D project variables
Fl_Button * photonLeft;
Fl_Button * photonRight;
Fl_Output * photonOutput;
Fl_Input *numPhotonHistories;
int simulationStatus = 0;
int simulationStatus3 = 0;
Fl_Button * photonLeft3;
Fl_Button * photonRight3;
Fl_Output * photonOutput3;

float scale = 0.1;
float manualLR = 0.0;
float glUp2D = 0.0;
float glLR2D = 0.0;

#include "utils.cxx"
#include "3d_opengl.cxx"
#include "2d_opengl.cxx"
#include "visual_win_3d_opengl.cxx"

/****************************************************************************************
 * Functions below are for
 * the operation of the GUI
 * *************************************************************************************/

//menu
Fl_Menu_Item menu_menu[] = {
 {"View/Edit Parameters", 0,  0, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0}
};

Fl_Menu_Bar * menuBar;
Fl_Menu_Button * menuButton;
Fl_Menu_Button * menuStart;

//subwindow
Fl_Window * subWindow;
Fl_Window * visualWindow;
Fl_Button * hideButton;
Fl_Button * saveButton;
Fl_Button * cancelButton;
int hideCounter = 0;

//inputs
Fl_Input *xRay;
Fl_Input *minD;
Fl_Input *maxD;
Fl_Input *bins;
Fl_Input *xDim;
Fl_Input *yDim;
Fl_Input *dThickness;
Fl_Input *colRad;
Fl_Input *colRefrac;
Fl_Input *interColRefrac;
Fl_Input *absorpFrac;
Fl_Input *bulkAbsorpCoeff;
Fl_Input *surRoughCoeff;
Fl_Input *minDistance;
Fl_Input *maxDistance;
Fl_Input *xLower;
Fl_Input *yLower;
Fl_Input *xUpper;
Fl_Input *yUpper;
Fl_Input *lightYield;
Fl_Input *pixelPitch;
Fl_Input *niSensor;
Fl_Input *flag;
Fl_Input *machine;

//button
Fl_Button * startButton;
Fl_Button * plus;
Fl_Button * minus;
Fl_Button * visual;

//output
Fl_Text_Display * output;
Fl_Text_Buffer *buff;

//window
Fl_Window * o;

//boxes
Fl_Box * b1;
Fl_Box * b2;
Fl_Box * b3;
Fl_Box * b4;
Fl_Box * b5;
Fl_Box * imageLabel;
Fl_Box * detectLabel;

//hybridMANTIS.out file reader
int fileCreated = 0;
FILE *outputfp;
char line [1000];

//myimage variables
int imageCBcount = 1;
int imageLRCBcount = 0;
Fl_PNG_Image * image;
Fl_Button * left;
Fl_Button * right;
Fl_Output * imageoutput;

//detect variables
Fl_PNG_Image * image2;
int detectCBcount = 1;
int detectLRCBcount = 0;
Fl_Button * detectLeft;
Fl_Button * detectRight;
Fl_Output * detectOutput;

//progress bar
Fl_Progress *progress;
int progressCounter = 0;

//temp
char temp[10000];
char path[1000];

#include "user_input.cxx"
#include "write_data.cxx"
#include "text_buffer.cxx"

/******************************
 * thread_task
 * 
 * Creates a new thread
 * to launch hybridMANTIS
 * and reads the output from
 * console.
 * ***************************/
void* thread_task(void* p) 
{
	FILE *console;
	
	console = popen("./visualMANTIS_ver1_0.x < penEasy_CsI_input.in", "r");
	while (fgets(path, 1000, console) != NULL)
	{
		append(path);
	}

	pclose(console);

     return 0;
}

/******************************
 * menu_cb
 * 
 * Shows the input subwindow
 * ***************************/
void menu_cb(void *)
{
	subWindow->show();
}

#include "my_detect_cb.cxx"
#include "my_image_cb.cxx"



/**********
 * visual_cb
 *********/
void visual_cb(Fl_Widget *) 
{
	printf("pressed\n");
	visualWindow->show();
}

/******************************
 * bcallback()
 * 
 * Callback corresponding to
 * "Start Simulation" button.
 * ***************************/
void bcallback(Fl_Widget *) 
{
	pid_t my_pid, parent_pid, child_pid;
	
	//if simulation hasn't started, set pressed = 1
	if(pressed == 0)
	{
		pressed = 1;
		
		//************************3d photon*******************//
		numHistStruct = atoi(numPhotonHistories->value());
		photonHistory = (struct histStruct *) malloc (numHistStruct * sizeof(struct histStruct));
		photonHistory3 = (struct histStruct *) malloc (10 * sizeof(struct histStruct));
	
		//write data files
		writeData();
		writePenEasy();

		FILE *console;	
		console = popen("rm *.dat *.png", "r");
		pclose(console);

		//run hybridMANTIS in separate thread
		Fl_Thread thread_id;
		fl_create_thread(thread_id, thread_task, (void *)NULL);
	}
}

/******************************
 * progress_cb()
 * 
 * Callback that reads the
 * deposition events data file
 * and updates the progress
 * bar.
 * ***************************/
void progress_cb(void *)
{
	FILE * depositionFP;
	int pastHist;
	char line[100];
	float value;
	
	depositionFP = fopen("CsI_deposition_events.dat", "r");
	
	//if simulation has started and deposition events exists
	if(depositionFP != NULL && pressed == 1)
	{
		//reads in number of hist simulated
		while (fgets(line, 100, depositionFP) != NULL)
		{
			fscanf(depositionFP, "#  Past hist =        %d ;", &pastHist);
		}
		fclose(depositionFP);
		value = pastHist / (atoi(xRay->value()) * 1.0);
		
		// update progress bar's label
		progress->value(value);
		char percent[100];
		sprintf(percent, "%.2f%% (%d / %s)", (pastHist/atof(xRay->value())) * 100.0, pastHist, xRay->value());
		progress->label(percent);              
		Fl::check();  
	}
	Fl::add_timeout(1.0, progress_cb);
}

/**************************************************************
 * initMainWindow()
 * 
 * Initializes and creates the main GUI window for hybridMANTIS.
 * ***********************************************************/
void initMainWindow()
{
	//initializes main window
	o = new Fl_Window(1740, 945, "visualMANTIS - v0.1");
    o->color(48);
    o->selection_color((Fl_Color)43);
	
	//Boxes for different graphical views	
		//detect
		b2 = new Fl_Box(1150, 55, 580, 445, "");
    	b2->box(FL_ENGRAVED_BOX);
      	b2->align(Fl_Align(FL_ALIGN_TOP|FL_ALIGN_INSIDE));

		//myimage
		b3 = new Fl_Box(1150, 500, 580, 445, "");
      	b3->box(FL_ENGRAVED_BOX);
      	b3->align(Fl_Align(FL_ALIGN_TOP|FL_ALIGN_INSIDE));
	
		//output dump
		b4 = new Fl_Box(0, 630, 1150, 3000, "output results");
      	b4->box(FL_ENGRAVED_BOX);
      	b4->labelfont(1);
      	b4->align(Fl_Align(FL_ALIGN_TOP|FL_ALIGN_INSIDE));
	
	//progress bar
		progress = new Fl_Progress(420,0,300,20);
		progress->minimum(0);                      // set progress range to be 0.0 ~ 1.0
		progress->maximum(1);
		progress->color(0x88888800);               // background color
		progress->selection_color(0x4444ff00);     // progress bar color
		progress->labelcolor(FL_WHITE);            // percent text color
    
	//buttons
		startButton = new Fl_Button(145, 0, 145, 20, "Start Simulation");
    	startButton->callback((Fl_Callback*) bcallback);
    
    //visual
		visual = new Fl_Button(290, 0, 125, 20, "Mini-Visualization");
		visual->callback((Fl_Callback*) visual_cb);
		
	//output
		output = new Fl_Text_Display(10, 660, 1130, 280, "");
		buff = new Fl_Text_Buffer();
		output->buffer(buff);
		buff->text("");
		
	//image controller
		left = new Fl_Button(1150+40, 910, 100, 20, "<--");
		right = new Fl_Button(1150+450, 910, 100, 20, "-->");
		imageoutput = new Fl_Output(1150+230, 910, 160, 20, "");		
		left->callback((Fl_Callback*) left_cb);
		right->callback((Fl_Callback*) right_cb);
	
	//detect controller
		detectLeft = new Fl_Button(1150+40, 460, 100, 20, "<--");
		detectRight = new Fl_Button(1150+450, 460, 100, 20, "-->");
		detectOutput = new Fl_Output(1150+230, 460, 160, 20, "");		
		detectLeft->callback((Fl_Callback*) detectLeft_cb);
		detectRight->callback((Fl_Callback*) detectRight_cb);
		
	//3d image controller
		photonLeft = new Fl_Button(620, 600, 100, 20, "<--");
		photonRight = new Fl_Button(1000, 600, 100, 20, "-->");
		photonOutput = new Fl_Output(800, 600, 160, 20, "");
		photonLeft->callback((Fl_Callback*) photonLeft_cb);
		photonRight->callback((Fl_Callback*) photonRight_cb);
		photonLeft->deactivate();
		photonRight->deactivate();
		
    //menu
		menuBar = new Fl_Menu_Bar(0, 0, 145, 20);
		menuButton = new Fl_Menu_Button(0, 0, 145, 20, "Options");
		menuButton->menu(menu_menu);
		menu_menu[0].callback((Fl_Callback*) menu_cb);
	
	//labels
		imageLabel = new Fl_Box(1150+200,475,200,100, "Point Response Function");
		imageLabel->labelsize(20);
		imageLabel->align(FL_ALIGN_CENTER);
		
		detectLabel = new Fl_Box(1150+200,25,200,100, "Pulse Height Spectrum");
		detectLabel->labelsize(20);
		detectLabel->align(FL_ALIGN_CENTER);
		
	o->resizable(o);	
	o->end();
	o->show();

	//1 second callbacks for background processing
	Fl::add_timeout(1.0, myimage_cb);
	Fl::add_timeout(1.0, detect_cb);
	Fl::add_timeout(1.0, progress_cb);
}

/**************************************************************
 * init3DprojectionWindow()
 * 
 * A sub OpenGL window is created to handle the 3D view of
 * photons.
 * ***********************************************************/
void init3DprojectionWindow()
{
	//o begin
	o->begin();
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	
	glutInitWindowPosition(620,60);
	glutInitWindowSize(500, 500);
	mainWindow = glutCreateWindow("Lighthouse3D - GLUT Tutorial");

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND); //Enable alpha blending
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function

	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutTimerFunc(1000, render, 0);

	// here are the two new functions;
	glutMotionFunc (glutMotion);
	glutMouseFunc(mouseScroll);		
	o->end();
}

/**************************************************************
 * init2DprojectionWindow()
 * 
 * A sub OpenGL window is created to handle the 2D view of
 * photons.
 * ***********************************************************/
void init2DprojectionWindow()
{
	//o begin
	o->begin();
	//glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(60,60);
	glutInitWindowSize(500, 500);
	mainWindow2 = glutCreateWindow("2D Projection");
	//glEnable(GL_DEPTH_TEST);
	//glEnable(GL_NORMALIZE);
	//glEnable(GL_COLOR_MATERIAL);
	//glEnable(GL_BLEND); //Enable alpha blending
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function
	
	// register callbacks
	glutDisplayFunc(renderScene2);
	glutReshapeFunc(changeSize2);
	glutMouseFunc(mouseScroll2);
	glutSpecialFunc(pressKey);
	//glutTimerFunc(1120, render2, 0);


	o->end();
}

/**************************************************************
 * initVisualWindow()
 * ***********************************************************/
void initVisualWindow()
{
	visualWindow = new Fl_Window(0,0, 420, 450, "Mini-Visual Window");
	visualWindow->color(48);
    visualWindow->selection_color((Fl_Color)43);
	visualWindow->begin();
	//3d image controller
	photonLeft3 = new Fl_Button(10, 420, 100, 20, "<--");
	photonRight3 = new Fl_Button(290, 420, 100, 20, "-->");
	photonOutput3 = new Fl_Output(120, 420, 160, 20, "");
	photonLeft3->callback((Fl_Callback*) photonLeft_cb3);
	photonRight3->callback((Fl_Callback*) photonRight_cb3);
	photonLeft3->deactivate();
	photonRight3->deactivate();
		
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
	
	glutInitWindowPosition(5, 5);
	glutInitWindowSize(400, 400);
	mainWindow3 = glutCreateWindow("");

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND); //Enable alpha blending
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function

	// register callbacks
	glutDisplayFunc(renderScene3);
	glutReshapeFunc(changeSize3);
	glutTimerFunc(1100, render3, 0);

	// here are the two new functions;
	glutMotionFunc (glutMotion3);
	glutMouseFunc(mouseScroll3);		
	
	visualWindow->resizable(visualWindow);
	visualWindow->end();
	visualWindow->hide();
}

/**************************************************************
 * main()
 * 
 * The main function that initializes and starts the GUI.
 * ***********************************************************/
int main(int argc, char **argv)
{
	initParameterWindow();
	initMainWindow();
	init3DprojectionWindow();
	init2DprojectionWindow();
	initVisualWindow();

	return Fl::run();
}
