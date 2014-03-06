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
//	Filename:	visual_win_3d_opengl.h
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file contains the code to handle the 3D photon visualization window.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//colors for the cylinders
float lineColors3[10][3] = 
{
	{1.0, 1.0, 1.0},
	{1.0, 1.0, 0.0},
	{1.0, 1.0, 0.5},
	{1.0, 1.0, 0.75},
	{1.0, 0.0, 0.25},
	{1.0, 0.0, 1.0},
	{1.0, 0.25, 1.0},
	{1.0, 0.5, 1.0},
	{1.0, 0.75, 1.0},
	{1.0, 0.9, 1.0}
};

float cylinderColors3[10][3] = 
{
	{0.0f, 0.0f, 1.0f},
	{0.0f, 1.0f, 0.0f},
	{1.0f, 0.0f, 0.0f},
	{0.0f, 1.0f, 1.0f},
	{1.0f, 1.0f, 0.0f},
	{1.0f, 1.0f, 1.0f},
	{0.2f, 0.0f, 0.0f},
	{0.4f, 0.5f, 0.3f},
	{0.6f, 0.0f, 0.4f},
	{0.8f, 0.8f, 0.0f}
};

// angle of rotation for the camera direction
float angle3= 0.0f;

// actual vector representing the camera's direction
float lx3=0.0f,lz3=-1.0f, ly3=0.0f;

// XZ position of the camera
float x3=1.0f, y3=1.0f, z3=1200.0f;

// the key states. These variables will be zero
//when no key is being presses
float deltaAngle3 = 0.0f;
float deltaMove3 = 0;
int xOrigin3 = -1;
int keyPress3 = 0;

/* old position of the mouse */
int oldX3 = -13;
int oldY3 = -13;

/* mouse state, UP or DOWN */
int mState3 = UP;

/* current axis of rotation */
int axisRot3 = X_AXIS;

/* global rotation, for use with the mouse */
vec3_t gRot3 = {0, 0, 0};

/******************************
 * changeSize
 * 
 * Initiates size of
 * 3D projection first.
 * ***************************/
void changeSize3(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;

	float ratio =  w * 1.0 / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(90.0f, ratio, 0.1f, 10000.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);

	glutSetWindow(mainWindow3);
	glutPostRedisplay();
}

/******************************
 * computePos
 * 
 * Computes position of camera
 * for  3D scene
 * ***************************/
void computePos3(float deltaMove) 
{
	x3 += deltaMove3 * lx3 * 0.1f;
	z3 += deltaMove3 * lz3 * 0.1f;

	glutSetWindow(mainWindow3);
	glutPostRedisplay();
}

/******************************
 * computeDir
 * 
 * Computes the direction
 * for the camera
 * ***************************/
void computeDir3(float deltaAngle) 
{
	angle3 += deltaAngle;
	lx3 = sin(angle3);
	lz3 = -cos(angle3);

	glutSetWindow(mainWindow3);
	glutPostRedisplay();
}

/******************************
 * drawLowerNet3
 * 
 * Draws the grid for the 3D projection
 * ***************************/
void drawLowerNet3(int size, int lx, int lz)
{
	int xc, zc;
	
	glColor3f(0.9f, 0.9f, 0.9f);
	glPushMatrix();
		glPushAttrib(GL_ENABLE_BIT); 
		glLineStipple(1, 0xAAAA);
		glEnable(GL_LINE_STIPPLE);
	
		//draws lines between vertices
		glBegin(GL_LINES);
		xc = 0;
		vertx = -size / 2.0 + xc / (GLfloat)(lx-1)*size;
		verty = 0.0;
		vertz = size / 2.0;
			
		glVertex3f(0.0 - atoi(def_xUpper) / 2.0, 0.0, 0.0 - atoi(def_xUpper) / 2.0);
		glVertex3f(size / 2.0 - atoi(def_xUpper) / 2.0, 0.0, 0 - atoi(def_xUpper) / 2.0);
			
		glVertex3f(0.0 - atoi(def_xUpper) / 2.0, 0.0, 0.0 - atoi(def_xUpper) / 2.0);
		glVertex3f(0.0 - atoi(def_xUpper) / 2.0, 0.0, size / 2.0 - atoi(def_xUpper) / 2.0);
			
		glVertex3f(size / 2.0 - atoi(def_xUpper) / 2.0, 0.0, 0 - atoi(def_xUpper) / 2.0);
		glVertex3f(size / 2.0 - atoi(def_xUpper) / 2.0, 0.0, size / 2.0 - atoi(def_xUpper) / 2.0);
			
		glVertex3f(0.0 - atoi(def_xUpper) / 2.0, 0.0, size / 2.0 - atoi(def_xUpper) / 2.0);
		glVertex3f(size / 2.0 - atoi(def_xUpper) / 2.0, 0.0, size / 2.0 - atoi(def_xUpper) / 2.0);
			
		glVertex3f(0.0 - atoi(def_xUpper) / 2.0, 0.0, size / 4.0 - atoi(def_xUpper) / 2.0);
		glVertex3f(size / 2.0 - atoi(def_xUpper) / 2.0, 0.0, size / 4.0 - atoi(def_xUpper) / 2.0);
			
		glVertex3f(size / 4.0 - atoi(def_xUpper) / 2.0, 0.0, 0.0 - atoi(def_xUpper) / 2.0);
		glVertex3f(size / 4.0 - atoi(def_xUpper) / 2.0, 0.0, size / 2.0 - atoi(def_xUpper) / 2.0);
		glEnd();
		glPopAttrib();
		
	glPopMatrix();
}

/******************************
 * photonLeftDynamic3
 * 
 * Gets called by left button
 * click to change 3D projection to
 * next data set. Displays the
 * data dynamically by the timer.
 * ***************************/
void photonLeftDynamic3(int filenum)
{
	int i, terminate;
	FILE *photonHistFP;
	char histname[100];
	float Xc, Yc, x, y, z;
	char phistLine [100];
	
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	if(photonHistFP != NULL && pressed == 1 && photonHistory3[glDataCounter3].dynamicCounter < photonHistory3[glDataCounter3].histCounter && filenum == glDataCounter3)
	{		
  		fgets(phistLine, 100, photonHistFP);
  		photonHistory3[glDataCounter3].histCounter = atoi(phistLine);
		photonHistory3[glDataCounter3].dynamicCounter += 1;
		
		for(i=0;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
		{
			fgets(phistLine, 100, photonHistFP);
			fscanf(photonHistFP, "%f %f %f %f %f %d", &x, &y, &z, &Xc, &Yc, &terminate);
			photonHistory3[glDataCounter3].x[i] = x;
			photonHistory3[glDataCounter3].y[i] = y;
			photonHistory3[glDataCounter3].z[i] = z+75.0;
			photonHistory3[glDataCounter3].Xc[i] = Xc;
			photonHistory3[glDataCounter3].Yc[i] = Yc;
			photonHistory3[glDataCounter3].terminated[i] = terminate;
			
			if((abs(x - Xc) > 6.0 || abs(y - Yc) > 6.0) && (abs(z) == 75.0))
			{
				photonHistory3[glDataCounter3].ceiling[i] = 1;
			}
			else
			{
				photonHistory3[glDataCounter3].ceiling[i] = 0;
			}
		}
  		fclose(photonHistFP);
  		
  		photonOutput3->value(histname);
		photonOutput3->redraw();
		
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
		
		glutTimerFunc(11, photonLeftDynamic3, filenum);
	}
}

/******************************
 * photonRightDynamic
 * 
 * Gets called by right button
 * click to change 3D projection to
 * next data set. Displays the
 * data dynamically by the timer.
 * ***************************/
void photonRightDynamic3(int filenum)
{
	int i, terminate;
	FILE *photonHistFP;
	char histname[100];
	float Xc, Yc, x, y, z;
	char phistLine [100];
	
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	if(photonHistFP != NULL && pressed == 1 && photonHistory3[glDataCounter3].dynamicCounter < photonHistory3[glDataCounter3].histCounter && filenum == glDataCounter3)
	{		
  		fgets(phistLine, 100, photonHistFP);
  		photonHistory3[glDataCounter3].histCounter = atoi(phistLine);
		photonHistory3[glDataCounter3].dynamicCounter += 1;
		
		for(i=0;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
		{
			fgets(phistLine, 100, photonHistFP);
			fscanf(photonHistFP, "%f %f %f %f %f %d", &x, &y, &z, &Xc, &Yc, &terminate);
			photonHistory3[glDataCounter3].x[i] = x;
			photonHistory3[glDataCounter3].y[i] = y;
			photonHistory3[glDataCounter3].z[i] = z+75.0;
			photonHistory3[glDataCounter3].Xc[i] = Xc;
			photonHistory3[glDataCounter3].Yc[i] = Yc;
			photonHistory3[glDataCounter3].terminated[i] = terminate;
			
			if((abs(x - Xc) > 6.0 || abs(y - Yc) > 6.0) && (abs(z) == 75.0))
			{
				photonHistory3[glDataCounter3].ceiling[i] = 1;
			}
			else
			{
				photonHistory3[glDataCounter3].ceiling[i] = 0;
			}
		}
  		fclose(photonHistFP);
  		
  		photonOutput3->value(histname);
		photonOutput3->redraw();
		
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
		glutTimerFunc(11, photonRightDynamic3, filenum);
	}
}

/******************************
 * rewind3 to 0
 * 
 * Callback that handles
 * the left button click to
 * change the 3D photon history.
 * ***************************/
void rewind3(int num)
{
	FILE *photonHistFP;
	char histname[100];
	
	glDataCounter3 = 0;
	
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	if(photonHistFP != NULL && pressed == 1)
	{	
		printf("photonLeft: %d\n", glDataCounter3);
		fclose(photonHistFP);

		photonHistory3[glDataCounter3].dynamicCounter = 0;
		glutTimerFunc(11, photonLeftDynamic3, glDataCounter3);
		
		photonOutput3->value(histname);
		photonOutput3->redraw();
		photonLeft3->activate();
		photonRight3->activate();
	}
}

/******************************
 * render
 * 
 * Reads in the photon
 * data from file, and calls
 * renderScene functionto display
 * ***************************/
void render3(int num)
{
	int i, terminate;
	FILE *photonHistFP;
	char histname[100];
	float Xc, Yc, x, y, z;
	char phistLine [100];
		
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	//printf("glDataCounter: %d histname: %s\n", glDataCounter, histname);
	
	if(photonHistFP != NULL && pressed == 1)
	{		
  		fgets(phistLine, 100, photonHistFP);
  		photonHistory3[glDataCounter3].histCounter = atoi(phistLine);
  		photonHistory3[glDataCounter3].numTerminates = 0;
  				
		for(i=0;i<photonHistory3[glDataCounter3].histCounter;i++)
		{
			fgets(phistLine, 100, photonHistFP);
			fscanf(photonHistFP, "%f %f %f %f %f %d", &x, &y, &z, &Xc, &Yc, &terminate);
			photonHistory3[glDataCounter3].x[i] = x;
			photonHistory3[glDataCounter3].y[i] = y;
			photonHistory3[glDataCounter3].z[i] = z+75.0;
			photonHistory3[glDataCounter3].Xc[i] = Xc;
			photonHistory3[glDataCounter3].Yc[i] = Yc;
			photonHistory3[glDataCounter3].terminated[i] = terminate;
			
			if((abs(x - Xc) > 6.0 || abs(y - Yc) > 6.0) && (abs(z) == 75.0))
			{
				photonHistory3[glDataCounter3].ceiling[i] = 1;
			}
			else
			{
				photonHistory3[glDataCounter3].ceiling[i] = 0;
			}
			
			if(terminate > 0)
			{
				photonHistory3[glDataCounter3].numTerminates += 1;
			}
			
			if(i == (photonHistory3[glDataCounter3].histCounter-1) && photonHistory3[glDataCounter3].terminated[i] == 0)
			{
				photonHistory3[glDataCounter3].terminated[i] = 1;
				photonHistory3[glDataCounter3].numTerminates += 1;
			}
		}
  		fclose(photonHistFP);
  		
  		//printf("glDataCounter: %d %s histCounter: %d ", glDataCounter, histname, photonHistory3[glDataCounter].histCounter);

  		photonOutput3->value(histname);
		photonOutput3->redraw();
		
		glDataCounter3++;
		
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
		//d->redraw();
	}
	
	if(glDataCounter3 == 10)
	{
		//glDataCounter -= 1;
		//simulationStatus = 1;
		//glutSetWindow(mainWindow);
		//glutPostRedisplay();
	}
	else
	{
		glutTimerFunc(1100, render3, 0);
	}
	//else
	//{
		//glDataCounter = atoi(numPhotonHistories->value()) - 1;
		//simulationStatus = 1; //simulation over
		//glutTimerFunc(1000, render, 0); 
	//}
}

/******************************
 * renderScene
 * 
 * renders the 3D projection
 * ***************************/
void renderScene3(void) 
{	
	if (deltaMove3)
	{
		computePos3(deltaMove3);
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
	}

	if(deltaAngle3 && xOrigin3 < 0 && keyPress3 == 1)
	{
		computeDir3(deltaAngle3);
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
	}
	
	glMatrixMode(GL_MODELVIEW);

	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();

	// Set the camera
	gluLookAt(x3, y3, z3,
			x3+lx3, 1.0f, z3+lz3,
			0.0f, 1.0f, 0.0f);
    	
    	
    glRotatef (gRot3[0], 1.0, 0.0, 0.0);
    glRotatef (gRot3[1], 0.0, 1.0, 0.0);		
	
	//1.000000 1.000000 1248.400269 1.000000 1.000000 1247.400269 0.000000 1.000000 0.000000 91.800049 -0.000005
	//printf("gluLookAt: x:%f y:%f z:%f x+lx:%f %f z+lz:%f %f %f %f gRot[0]:%f gRot[1]:%f\n", x, y, z, x+lx, 1.0f, z+lz, 0.0f, 1.0f, 0.0f, gRot[0], gRot[1]);
	
	GLfloat size = atoi(def_xUpper) * 2;
	
	//drawGround();
	//drawAxis();
	//drawCeil();
	
	char temp[100];
	
	//lower net
	drawLowerNet3(size, 5, 5);
	glColor3f(1.0, 1.0, 0.0);
	sprintf(temp, "%d", 0);
	renderBitmapString(-atoi(def_xUpper)/2, verty, -1 * atoi(def_xUpper)/2, (void *)GLUT_BITMAP_TIMES_ROMAN_10, temp);
	sprintf(temp, "%d", atoi(def_xUpper)/2);
	renderBitmapString(0, verty, -atoi(def_xUpper)/2, (void *)GLUT_BITMAP_TIMES_ROMAN_10, temp);
	sprintf(temp, "%d", atoi(def_xUpper));
	renderBitmapString(atoi(def_xUpper)/2, verty, -atoi(def_xUpper)/2, (void *)GLUT_BITMAP_TIMES_ROMAN_10, temp);
	glColor3f(1.0, 0.0, 0.0);
	renderBitmapString(0, verty, -atoi(def_xUpper)/2 - 100, (void *) GLUT_BITMAP_TIMES_ROMAN_10, (char *)"x ( um )");

	glColor3f(1.0, 1.0, 0.0);	
	sprintf(temp, "%d", atoi(def_xUpper)/2);
	renderBitmapString(-atoi(def_xUpper)/2, verty, 0, (void *)GLUT_BITMAP_TIMES_ROMAN_10, temp);
	sprintf(temp, "%d", atoi(def_xUpper));
	renderBitmapString(-atoi(def_xUpper)/2, verty, atoi(def_xUpper)/2, (void *)GLUT_BITMAP_TIMES_ROMAN_10, temp);
	glColor3f(1.0, 0.0, 0.0);
	renderBitmapString(-atoi(def_xUpper)/2 - 150, verty, 0, (void *)GLUT_BITMAP_TIMES_ROMAN_10, (char *)"y ( um )");	

	float ax, ay, bx, by, da, db;
	//da, db are angles from 0,0
	da = 90.0;
	db = 90.0;

	//parametric function for circumference, 2.0 is radius
	ax = 0.0 + 2.0 * Cos(da);
	ay = 0.0 + 2.0 * Sin(da);
	bx = 10.0 + 2.0 * Cos(db);
	by = 0.0 + 2.0 * Sin(db);

	int ele = glDataCounter3 - 1;
	
	if(simulationStatus3 == 0 && pressed == 1)
	{
		int i, j, newPhoton;
		
		newPhoton = 0;
		i = 0;
		for(j=0;j<photonHistory3[ele].numTerminates; j++)
		{
			glColor3f(lineColors[j][0], lineColors[j][1], lineColors[j][2]);
			//glColor3f(1.0f, 1.0f, 0.0f);
			glBegin(GL_LINE_STRIP);//start drawing a line loop
			if(glDataCounter3 >= 0 && photonHistory3[ele].histCounter > 0)
			{
				for(i=newPhoton;i<photonHistory3[ele].histCounter;i++)
				{
					if(photonHistory3[ele].z[i] >= 0 && photonHistory3[ele].z[i] <= 150.0)
					{
						glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2), photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2));
					}
					if(photonHistory3[ele].terminated[i] > 0)
					{
						newPhoton = i;
						break;
					}
				}
			}
			glEnd();//end drawing of line loop
		}
	
		glColor3f(1.0f, 1.0f, 0.0f);
		glPointSize(5.0f);
		glBegin(GL_POINTS); //starts drawing of points
		if(glDataCounter3 >= 0 && photonHistory3[ele].histCounter > 0)
		{
			for(i=0;i<photonHistory3[ele].histCounter;i++)
			{
				if(photonHistory3[ele].z[i] >= 0 && photonHistory3[ele].z[i] <= 150.0)
				{
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2), photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2));
				}
			}
		}
		glEnd();//end drawing of points

		newPhoton = 0;
		i = 0;
		for(j=0;j<photonHistory3[ele].numTerminates; j++)
		{
			if(glDataCounter3 >= 0 && photonHistory3[ele].histCounter > 0)
			{
				for(i=newPhoton;i<photonHistory3[ele].histCounter;i++)
				{
					glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
					glEnable(GL_DEPTH_TEST);
					glEnable(GL_NORMALIZE);
					glEnable(GL_COLOR_MATERIAL);
					glEnable(GL_BLEND); //Enable alpha blending
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function
					glPushMatrix();
					//printf("%d %f %f\n", i, photonHistory3[ele-1].Xc[i], photonHistory3[ele-1].Yc[i]);
					glTranslatef(photonHistory3[ele].Xc[i] - (atoi(def_xUpper)/ 2), 0.0, photonHistory3[ele].Yc[i] - (atoi(def_xUpper)/ 2));
					glPushMatrix();
					glDisable(GL_TEXTURE_2D);
					GLUquadricObj * qobj = gluNewQuadric();
		
					//rotates cylinder to be upright
					glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
		
					//draws bottom disk
					glColor4f(cylinderColors[j][0], cylinderColors[j][1], cylinderColors[j][2], ALPHA);
					gluDisk(qobj, 0.0, atof(def_colRad), 32, 1);

					//draws actual cylinder
					glDisable(GL_TEXTURE_2D);
					glColor4f(cylinderColors[j][0], cylinderColors[j][1], cylinderColors[j][2], ALPHA);
					gluCylinder(qobj, atof(def_colRad), atof(def_colRad), 150.0f, 32, 1);

					//draws top disk
					glDisable(GL_TEXTURE_2D);
					glTranslatef(0.0f, 0.0f, 150.0f);
					glColor4f(1.0f, 0.0f, 0.0f, ALPHA);
					gluDisk(qobj, 0.0, atof(def_colRad), 32, 1);
					glPopMatrix();
					glPopMatrix();
				
					if(photonHistory3[ele].terminated[i] > 0)
					{
						newPhoton = i+1;
						break;
					}
				}
			}
		}

		if(glDataCounter3 >= 0 && photonHistory3[ele].histCounter > 0)
		{
			for(i=0;i<photonHistory3[ele].histCounter;i++)
			{
				if(photonHistory3[ele].z[i] == 0.00)
				{
					glDisable(GL_TEXTURE_2D);
					glColor4f(1.0f, 0.5f, 0.0f, ALPHA);
					glBegin(GL_POLYGON);//begin drawing of polygon
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glEnd();//end drawing of polygon
				}
				else if(photonHistory3[ele].z[i] == 150.00)
				{
					glDisable(GL_TEXTURE_2D);
					glColor4f(0.9f, 0.9f, 0.9f, ALPHA);
					glBegin(GL_POLYGON);//begin drawing of polygon
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glVertex3f(photonHistory3[ele].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[ele].z[i], photonHistory3[ele].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glEnd();//end drawing of polygon
				}
			}
		}
	}
	else if(simulationStatus3 == 1)
	{
		int i, j, newPhoton;
		
		newPhoton = 0;
		i = 0;
		for(j=0;j<photonHistory3[glDataCounter3].numTerminates; j++)
		{
			glColor3f(1.0, 1.0, 0.0);
			glBegin(GL_LINE_STRIP);//start drawing a line loop
			if(glDataCounter3 >= 0 && photonHistory3[glDataCounter3].histCounter > 0)
			{
				for(i=newPhoton;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
				{
					if(photonHistory3[glDataCounter3].z[i] >= 0 && photonHistory3[glDataCounter3].z[i] <= 150.0)
					{
						glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2), photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2));
					}
					
					if(photonHistory3[glDataCounter3].terminated[i] > 0)
					{
						newPhoton = i+1;
						break;
					}
				}
			}
			glEnd();//end drawing of line loop
		}
		
		glColor3f(1.0, 1.0, 0.0);
		glPointSize(5.0f);
		glBegin(GL_POINTS); //starts drawing of points
		if(glDataCounter3 >= 0 && photonHistory3[glDataCounter3].histCounter > 0)
		{
			for(i=0;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
			{
				if(photonHistory3[glDataCounter3].z[i] >= 0 && photonHistory3[glDataCounter3].z[i] <= 150.0)
				{
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2), photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2));
				}
			}
		}
		glEnd();//end drawing of points*/
			
		newPhoton = 0;
		i = 0;
		for(j=0;j<photonHistory3[glDataCounter3].numTerminates; j++)
		{
			if(glDataCounter3 >= 0 && photonHistory3[glDataCounter3].histCounter > 0)
			{
				for(i=newPhoton;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
				{
					//printf("%f %f\n", photonHistory3[glDataCounter3].Xc[i], photonHistory3[glDataCounter3].Yc[i]);
					glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
					glEnable(GL_DEPTH_TEST);
					glEnable(GL_NORMALIZE);
					glEnable(GL_COLOR_MATERIAL);
					glEnable(GL_BLEND); //Enable alpha blending
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function
					glPushMatrix();
					glTranslatef(photonHistory3[glDataCounter3].Xc[i] - (atoi(def_xUpper)/ 2), 0.0, photonHistory3[glDataCounter3].Yc[i] - (atoi(def_xUpper)/ 2));
					glPushMatrix();
					glDisable(GL_TEXTURE_2D);
					GLUquadricObj * qobj = gluNewQuadric();
		
					//rotates cylinder to be upright
					glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
		
					//draws bottom disk
					glColor4f(0.0f, 0.0f, 1.0f, ALPHA);
					gluDisk(qobj, 0.0, 5.1, 32, 1);
	
					//draws actual cylinder
					glDisable(GL_TEXTURE_2D);
					glColor4f(cylinderColors[j][0], cylinderColors[j][1], cylinderColors[j][2], ALPHA);
					gluCylinder(qobj, 5.1f, 5.1f, 150.0f, 32, 1);

					//draws top disk
					glDisable(GL_TEXTURE_2D);
					glTranslatef(0.0f, 0.0f, 150.0f);
					glColor4f(1.0f, 0.0f, 0.0f, ALPHA);
					gluDisk(qobj, 0.0, 5.1, 32, 1);
					glPopMatrix();
					glPopMatrix();
				
					if(photonHistory3[glDataCounter3].terminated[i] > 0)
					{
						newPhoton = i+1;
						break;
					}
				}
			}
		}
		
		
		if(glDataCounter3 >= 0 && photonHistory3[glDataCounter3].histCounter > 0)
		{
			for(i=0;i<photonHistory3[glDataCounter3].dynamicCounter;i++)
			{
				if(photonHistory3[glDataCounter3].z[i] == 0.00)
				{
					glDisable(GL_TEXTURE_2D);
					glColor4f(1.0f, 0.5f, 0.0f, ALPHA);
					glBegin(GL_POLYGON);//begin drawing of polygon
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glEnd();//end drawing of polygon
				}
				else if(photonHistory3[glDataCounter3].z[i] == 150.00)
				{
					glDisable(GL_TEXTURE_2D);
					glColor4f(0.9f, 0.9f, 0.9f, ALPHA);
					glBegin(GL_POLYGON);//begin drawing of polygon
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) - 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) + 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glVertex3f(photonHistory3[glDataCounter3].x[i] - (atoi(def_xUpper)/ 2) - 15.0, photonHistory3[glDataCounter3].z[i], photonHistory3[glDataCounter3].y[i] - (atoi(def_xUpper)/ 2) + 15.0);//first vertex
					glEnd();//end drawing of polygon
				}
			}
		}
	}

	glutSwapBuffers();

	if(simulationStatus3 == 0 && pressed == 1 && glDataCounter3 == 10)
	{
		rewind3(0);
		simulationStatus3 = 1;
	}
} 

/******************************
 * glutMotion
 * 
 * Calculates the mouse
 * movement.
 * ***************************/
void glutMotion3(int x, int y) 
{
	if (simulationStatus3 == 1)
	{
		if (mState3 == DOWN) 
		{
			gRot3[0] -= ((oldY3 - y) * 180.0f) / 100.0f;
			gRot3[1] -= ((oldX3 - x) * 180.0f) / 100.0f;
			//clamp (gRot);
			glutPostRedisplay ();
			//printf("%d %d %d %d\n", oldX, oldY, x, y);
		}	 
		oldX3 = x; 
		oldY3 = y;
		glutSetWindow(mainWindow3);
		glutPostRedisplay();
	}
}

/******************************
 * mouseScroll
 * 
 * Handles zoom function for
 * 3D projection
 * ***************************/
void mouseScroll3(int button, int state, int x, int y)
{
	if (simulationStatus3 == 1)
	{
		// Wheel reports as button 3(scroll up) and button 4(scroll down)
		if(x3 >= 0 && x3 <= 500 && y3 >= 0 && y3 <= 500)
		{
			if (button == 3) // It's a wheel event
			{	
				deltaMove3 = 88.0f;
			}
			else if(button == 4)
			{
				deltaMove3 = -88.0f;
			}
			else
			{
				deltaMove3 = 0.0f;
				// only start motion if the left button is pressed
				if (button == GLUT_LEFT_BUTTON) 
				{
					// when the button is released
					if (state == GLUT_UP) 
					{
						angle3 += deltaAngle3;
						xOrigin3 = -1;
					}
					else  
					{// state = GLUT_DOWN
						xOrigin3 = x;
						mState3 = DOWN;
						oldX3 = x;
						oldY3 = y;
					}
				}
			}
			
			glutSetWindow(mainWindow3);
			glutPostRedisplay();
		}
	}
}

/******************************
 * photonLeft_cb
 * 
 * Callback that handles
 * the left button click to
 * change the 3D photon history.
 * ***************************/
void photonLeft_cb3(void *)
{
	FILE *photonHistFP;
	char histname[100];
	
	glDataCounter3 --;
	
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	if(photonHistFP != NULL && pressed == 1)
	{	
		printf("photonLeft: %d\n", glDataCounter3);
		fclose(photonHistFP);

		photonHistory3[glDataCounter3].dynamicCounter = 0;
		glutTimerFunc(11, photonLeftDynamic3, glDataCounter3);
		
		photonOutput3->value(histname);
		photonOutput3->redraw();
	}
	else
	{
		glDataCounter3 ++;
	}
}

/********************************
 * photonRight_cb
 * 
 * Callback that handles
 * the right button click to
 * change the 3D photon history.
 * ******************************/
void photonRight_cb3(void *)
{
	FILE *photonHistFP;
	char histname[100];
	
	glDataCounter3 ++;
	
	sprintf(histname, "photon_hist_1_%d.dat", glDataCounter3);
	photonHistFP = fopen(histname, "r");
	
	//if the data file exists and simulation started
	if(photonHistFP != NULL && pressed == 1)
	{	
		printf("photonRight: %d\n", glDataCounter3);
		fclose(photonHistFP);
		photonHistory3[glDataCounter3].dynamicCounter = 0;
		glutTimerFunc(11, photonRightDynamic3, glDataCounter3);
		
		photonOutput3->value(histname);
		photonOutput3->redraw();
	}
	else
	{
		glDataCounter3 --;
	}
}
