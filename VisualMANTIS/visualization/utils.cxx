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
//	Filename:	utils.cxx
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file has utility files that are used throughout the visualization.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************
 * drawCylinder
 * 
 * Draws cylinder on the 3D
 * scene
 * ***************************/
void drawCylinder() 
{
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
		glColor4f(0.0f, 0.0f, 1.0f, ALPHA);
	 	gluCylinder(qobj, 5.1f, 5.1f, 150.0f, 32, 1);

		//draws top disk
		glDisable(GL_TEXTURE_2D);
		glTranslatef(0.0f, 0.0f, 150.0f);
		glColor4f(1.0f, 0.0f, 0.0f, ALPHA);
	 	gluDisk(qobj, 0.0, 5.1, 32, 1);
	glPopMatrix();
}



/******************************
 * 
 * 
 * ***************************/
void drawGround()
{
	glColor3f(0.9f, 0.9f, 0.9f);
	glBegin(GL_QUADS);
		glVertex3f(-30.0f, 0.0f, -30.0f);
		glVertex3f(-30.0f, 0.0f,  30.0f);
		glVertex3f( 30.0f, 0.0f,  30.0f);
		glVertex3f( 30.0f, 0.0f, -30.0f);
	glEnd();
}

/******************************
 * 
 * 
 * ***************************/
void drawCeil(int x, int y, int z)
{
	glDisable(GL_TEXTURE_2D);
	glColor4f(0.9f, 0.9f, 0.9f, ALPHA);
	glBegin(GL_POLYGON);//begin drawing of polygon
      glVertex3f(-1000.0f,150.0f,0.0f);//first vertex
      glVertex3f(1000.0f,150.0f,0.0f);//second vertex
      glVertex3f(1000.0f,150.0f,1000.0f);//third vertex
      glVertex3f(-1000.0f,150.0f,1000.0f);//fourth vertex
    glEnd();//end drawing of polygon
}

/******************************
 * drawAxis
 * 
 * not in use yet
 * ***************************/
void drawAxis (void)
{
    glColor3f (0.5, 0.5, 0.5);
    glBegin (GL_LINES);
        glColor3f (0.5, 0.0, 0.0);
        glVertex3f (-500.0, 0.0, 0.0);
        glVertex3f (500.0, 0.0, 0.0);

        glColor3f (0.0, 0.5, 0.0);
        glVertex3f (0.0, 500.0, 0.0);
        glVertex3f (0.0, -500.0, 0.0);

        glColor3f (0.0, 0.0, 0.5);
        glVertex3f (0.0, 0.0, -500.0);
        glVertex3f (0.0, 0.0, 500.0);
    glEnd ();
}

/******************************
 * drawNet
 * 
 * not in use yet
 * ***************************/
void drawNet(GLfloat size, GLint LinesX, GLint LinesZ)
{
	int xc, zc;
	
	glPushAttrib(GL_ENABLE_BIT); 
	glLineStipple(1, 0xAAAA);
	glEnable(GL_LINE_STIPPLE);

	glBegin(GL_LINES);
	for (xc = 0; xc < LinesX; xc++)
	{
		glVertex3f(	-size / 2.0 + xc / (GLfloat)(LinesX-1)*size, 0.0, size / 2.0);
		glVertex3f(	-size / 2.0 + xc / (GLfloat)(LinesX-1)*size, 0.0, size / -2.0);
	}
	for (zc = 0; zc < LinesZ; zc++)
	{
		glVertex3f(	size / 2.0, 0.0, -size / 2.0 + zc / (GLfloat)(LinesZ-1)*size);
		glVertex3f(	size / -2.0, 0.0, -size / 2.0 + zc / (GLfloat)(LinesZ-1)*size);
	}
	glEnd();
	glPopAttrib();
}



/******************************
 * drawLeftNet
 * 
 * not in use yet.
 * ***************************/
void drawLeftNet(int size, int lx, int lz)
{
	int xc, zc;
	
	glColor3f(1.0f, 0.0f, 0.0f);
	glPushMatrix();
		glRotatef(90.0f, 0.0f, 0.0f, 1.0f);
		glTranslatef(size/2, size/2, 0.0f);
		glPushAttrib(GL_ENABLE_BIT); 
		glLineStipple(1, 0xAAAA);
		glEnable(GL_LINE_STIPPLE);
		glBegin(GL_LINES);
		for (xc = 0; xc < lx; xc++)
		{
			glVertex3f(	-size / 2.0 + xc / (GLfloat)(lx-1)*size, 0.0, size / 2.0);
			glVertex3f(	-size / 2.0 + xc / (GLfloat)(lx-1)*size, 0.0, size / -2.0);
		}
		for (zc = 0; zc < lz; zc++)
		{
			glVertex3f(	size / 2.0, 0.0, -size / 2.0 + zc / (GLfloat)(lz-1)*size);
			glVertex3f(	size / -2.0, 0.0, -size / 2.0 + zc / (GLfloat)(lz-1)*size);
			
			if(zc == 0)
			{
				vertx = size/2.0;
				verty = 0.0;
				vertz = -size / 2.0 + zc / (GLfloat)(lz-1)*size;
			}
		}
		glEnd();
		glPopAttrib();
	glPopMatrix();
}

/******************************
 * drawRightNet
 * 
 * Not in use yet.
 * ***************************/
void drawRightNet(int size, int lx, int lz)
{
	int xc, zc;
	
	glColor3f(1.0, 1.0, 1.0);
	glPushMatrix();
		//glTranslatef(0.0f, -size/2, -size/2);
		//glPushAttrib(GL_ENABLE_BIT); 
		//glLineStipple(1, 0xAAAA);
		//glEnable(GL_LINE_STIPPLE);
	
		/*glBegin(GL_LINES);
		for (xc = 0; xc < lx; xc++)
		{
			glVertex3f(	-size / 2.0 + xc / (GLfloat)(lx-1)*size, 0.0, size / 2.0);
			glVertex3f(	-size / 2.0 + xc / (GLfloat)(lx-1)*size, 0.0, size / -2.0);
		}
		for (zc = 0; zc < lz; zc++)
		{
			glVertex3f(	size / 2.0, 0.0, -size / 2.0 + zc / (GLfloat)(lz-1)*size);
			glVertex3f(	size / -2.0, 0.0, -size / 2.0 + zc / (GLfloat)(lz-1)*size);

			if(zc == 0)
			{
				vertx = size/2.0;
				verty = 0.0;
				vertz = -size / 2.0 + zc / (GLfloat)(lz-1)*size;
			}
		}
		glEnd();*/
		glRotatef(90.0, 1.0f, 0.0, 0.0);		
		glRotatef(90.0, 1.0f, 0.0, 0.0);
		glDisable(GL_TEXTURE_2D);
		GLUquadricObj * qobj = gluNewQuadric();
		
		glTranslatef(20.0, -500.0, 0.0);
		
		//draws bottom disk
		glColor4f(1.0f, 1.0f, 1.0f, ALPHA);
		gluDisk(qobj, 0.0, 50.1, 32, 1);
	
		//draws actual cylinder
		/*glDisable(GL_TEXTURE_2D);
		glColor4f(0.0f, 0.0f, 1.0f, ALPHA);
	 	gluCylinder(qobj, 5.1f, 5.1f, 150.0f, 32, 1);

		//draws top disk
		glDisable(GL_TEXTURE_2D);
		glTranslatef(0.0f, 0.0f, 150.0f);
		glColor4f(1.0f, 0.0f, 0.0f, ALPHA);
	 	gluDisk(qobj, 0.0, 5.1, 32, 1);*/
	glPopMatrix();
	
}

/******************************
 * setOrthographicProjection
 * 
 * Initializes the projection
 * matrix for the scenes
 * ***************************/
void setOrthographicProjection() 
{
	// switch to projection mode
	glMatrixMode(GL_PROJECTION);
	// save previous matrix which contains the 
	//settings for the perspective projection
	glPushMatrix();
	// reset matrix
	glLoadIdentity();
	gluOrtho2D(0, 1024, 800,0);
	// set a 2D orthographic projection
	//gluOrtho2D(0, w, 0, h);
	//// invert the y axis, down is positive
	//glScalef(1, -1, 1);
	//// mover the origin from the bottom left corner
	//// to the upper left corner
	//glTranslatef(0, -h, 0);
	glMatrixMode(GL_MODELVIEW);
}

/******************************
 * resetPerspectiveProjection
 * 
 * ***************************/
void resetPerspectiveProjection() {
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

/******************************
 * renderBitmapString
 * 
 * Renders a string on the 
 * GL projections
 * ***************************/
void renderBitmapString(float x, float y, float z, void *font, char *string)
{
  
  char *c;
  glRasterPos3f(x, y, z);
  for (c=string; *c != '\0'; c++) {
    glutBitmapCharacter(font, *c);
  }
}

/******************************
 * renderSpacedBitmapString
 * 
 * not used yet
 * ***************************/
void renderSpacedBitmapString(float x, float y, float z, int spacing, void *font, char *string) {
  char *c;
  int x1=x;
  for (c=string; *c != '\0'; c++) {
	glRasterPos3f(x1,y,z);
    glutBitmapCharacter(font, *c);
	x1 = x1 + glutBitmapWidth(font,*c) + spacing;
  }
}

/******************************
 * renderVerticalBitmapString
 *
 * not used yet 
 * ***************************/
void renderVerticalBitmapString(float x, float y, float z, int bitmapHeight, void *font, char *string)
{
  
  char *c;
  int i;
  for (c=string,i=0; *c != '\0'; i++,c++) {
	glRasterPos3f(x, bitmapHeight-(i*2), z);
    glutBitmapCharacter(font, *c);
  }
}
