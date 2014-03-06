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
//	Filename:	2d_opengl.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file handles the 2-Dimensional top-down view of the optical photons. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// angle of rotation for the camera direction
float angle2 = 0.0f;

// actual vector representing the camera's direction
float lx2=0.0f,lz2=-1.0f, ly2=0.0f;

// XZ position of the camera
float x2=1.0f, y2=1.0f, z2=1200.0f;

// the key states. These variables will be zero
//when no key is being presses
float deltaAngle2 = 0.0f;
float deltaMove2 = 0;
int xOrigin2 = -1;
int keyPress2 = 0;

/* old position of the mouse */
int oldX2 = -13;
int oldY2 = -13;

/* mouse state, UP or DOWN */
int mState2 = UP;

/* current axis of rotation */
int axisRot2 = X_AXIS;

/******************************
 * changeSize2
 * 
 * Initiates size of
 * 2D projection first.
 * ***************************/
void changeSize2(int w, int h) {

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

	glutSetWindow(mainWindow2);
	glutPostRedisplay();
}

/******************************
 * computePos2
 * 
 * Computes position of camera
 * for 2D scene
 * ***************************/
void computePos2(float deltaMove2) 
{
	x2 += deltaMove2 * lx2 * 0.1f;
	z2 += deltaMove2 * lz2 * 0.1f;
	
	glutSetWindow(mainWindow2);
	glutPostRedisplay();
}

/******************************
 * drawLowerNet2
 * 
 * Draws the grid for the 2D projection
 * ***************************/
void drawLowerNet2(int size, int lx, int lz)
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
 * renderScene2
 * 
 * Renders the 2D projection
 * ***************************/
void renderScene2(void) 
{	
	if (deltaMove2)
	{
		computePos2(deltaMove2);
		glutSetWindow(mainWindow2);
		glutPostRedisplay();
	}
	
	glMatrixMode(GL_MODELVIEW);

	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();
	
	glTranslatef(glLR2D, glUp2D, 0.0);	
	gluLookAt(1.000000, 1.000000, z2,
			1.000000, 1.000000, z2+lz2,
			0.0f, 1.0f, 0.0f);
	glRotatef (90.0, 1.0, 0.0, 0.0);
    glRotatef (0.0, 0.0, 1.0, 0.0);
	
	GLfloat size = atoi(def_xUpper) * 2;
	
	char temp[100];
	
	//lower net
	drawLowerNet2(size, 5, 5);
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

	int i, j, counter, newPhoton;
		
	newPhoton = 0;
	i = 0;
	//int counter = glDataCounter-1;
	//for(counter=0;counter<glDataCounter2;counter++)
	//{	
	if(pressed == 1 && glDataCounter >= 0 && simulationStatus == 1)
	{
		counter = glDataCounter;
		for(i=newPhoton;i<photonHistory[counter].histCounter;i++)
		{	
			glPushMatrix();
			glTranslatef(photonHistory[counter].Xc[i] - (atoi(def_xUpper)/ 2), 0.0, photonHistory[counter].Yc[i]- (atoi(def_xUpper)/ 2));
				glPushMatrix();
				glDisable(GL_TEXTURE_2D);
							
				//rotates cylinder to be upright
				glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
		
				//draws bottom disk
				glColor4f(1.0f, 1.0f, 1.0f, ALPHA);
				float x,y;
				float radius = 5.1f;
				glBegin(GL_LINES);     
					x = (float)radius * cos(359 * PI/180.0f);
					y = (float)radius * sin(359 * PI/180.0f);
					for(int j = 0; j < 360; j++)
					{
						glVertex2f(x,y);
						x = (float)radius * cos(j * PI/180.0f);
						y = (float)radius * sin(j * PI/180.0f);
						glVertex2f(x,y);
					}
				glEnd();
					
				glPopMatrix();
			glPopMatrix();
		}
	}
	glutSwapBuffers();
} 

/******************************
 * pressKey
 * 
 * Handles up, down, left, 
 * right arrow buttons for
 * controlling 2D projection
 * ***************************/
void pressKey(int key, int xx, int yy) 
{	
	switch (key)
	{
		case GLUT_KEY_UP:
			glUp2D -= 100.0;
			break;
		case GLUT_KEY_DOWN: 
			glUp2D += 100.0;
			break;
		case GLUT_KEY_LEFT: 
			glLR2D += 100.0;
			break;
		case GLUT_KEY_RIGHT: 
			glLR2D -= 100.0;
			break;
	}
	deltaMove2 = 0.0;
	glutSetWindow(mainWindow2);
	glutPostRedisplay();
} 

/******************************
 * mouseScroll2
 * 
 * Handles zoom function
 * for 2D projection
 * ***************************/
void mouseScroll2(int button, int state, int x, int y)
{
	// Wheel reports as button 3(scroll up) and button 4(scroll down)
	if(x >= 0 && x <= 500 && y >= 0 && y <= 500)
	{
		if (button == 3) // It's a wheel event
		{	
			deltaMove2 = 388.0f;
		}
		else if(button == 4)
		{
			deltaMove2 = -388.0f;
		}
		else
		{
			deltaMove2 = 0.0f;
			// only start motion if the left button is pressed
			if (button == GLUT_LEFT_BUTTON) 
			{
				// when the button is released
				if (state == GLUT_UP) 
				{
					angle2 += deltaAngle2;
					xOrigin2 = -1;
				}
				else  
				{// state = GLUT_DOWN
					xOrigin2 = x;
					mState2 = DOWN;
					oldX2 = x;
					oldY2 = y;
				}
			}
		}
		glutSetWindow(mainWindow2);
		glutPostRedisplay();
	}
}
