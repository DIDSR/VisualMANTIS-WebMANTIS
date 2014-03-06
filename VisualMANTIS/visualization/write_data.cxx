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
//	Filename:	write_data.h
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This file writes data to text files for various uses throughout the visualization.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/******************************
 * writePenEasy
 * 
 * Writes to peneasy input
 * with user inputs
 * ***************************/
void writePenEasy()
{
	FILE *peasy;
	
	peasy = fopen("penEasy_CsI_input.in", "w");
	
	const char * text1 = 
	"# >>>> INPUT FILE FOR penEasy >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
#\n\
# CASE DESCRIPTION:\n\
#   penEasy with No EM fields!\n\
#\n\
# LAST UPDATE:\n\
#   2008-04-19 by JS\n\
#   2011-06-06 by DS\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR SECTION CONFIG\n\
#\n\
# * The simulation will be finished when ANY of the following conditions is\n\
#   fulfilled: i) the requested number of histories has been reached; ii) the\n\
#   alloted time has been exhausted; or iii) ALL the requested relative\n\
#   uncertainties have been reached.\n\
#\n\
# * Allotted time is interpreted as real time if it is a positive number.\n\
#   Otherwise, CPU (i.e. user) time is assumed.\n\
#\n\
# * Similarly, the update period is interpreted as real time if it is positive;\n\
#   otherwise it is the number of histories between updates. See the README file\n\
#   for a description of the tasks performed during each update. Note that the\n\
#   real time between updates cannot be larger than 50000 s.\n\
#\n\
# * Setting both random seeds equal to zero is an indication that they must be read\n\
#   from an external file. The name of the file containing the seeds is read from\n\
#   an additional input line that must be introduced after that of 'RANDOM SEEDS'.\n\
#   This new line must have the following format:\n\
#(rngseed.in                    ) RANDOM SEEDS FILE NAME (*** 30 characters ***)\n\
#   In this example, the file 'rngseed.in' is opened and read. The two seeds must\n\
#   be on the first line of rngseed.in, separated by one or more blanks. This\n\
#   feature is useful when using parallel computing (see README file).\n\
#\n\
[SECTION CONFIG v.2006-08-01]\n";
	
	const char * temp = xRay->value();
	int strl;
	int start;
	int i;
	char temp2[100];
	
	strl = strlen(temp);
	start = temp[0] - '0';
	
	//sprintf(temp2, "%1.2e", atof(xRay->value()));
	sprintf(temp2, "%d.%se%d", start, &temp[1], strl-1);
	//printf("%s %d %d %s\n", temp, strl, start, temp2);
	
	fprintf(peasy, "%s", text1);
	fprintf(peasy, "  %s                            NO. OF HISTORIES (<1e15)\n", temp2);

	const char * text2 =
"  5.0e15                           ALLOTTED TIME (s) (+ FOR REAL TIME; - FOR CPU TIME)\n\
  50000.0                          UPDATE INTERVAL (+ FOR REAL TIME (s) < 50000; - FOR HISTORIES)\n\
 3141592   20090625                INITIAL RANDOM SEEDS\n\
[END OF CONFIG SECTION]\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR SOURCE SECTIONS\n\
#\n\
# * Details on the features and configuration of each source are provided in their\n\
#   accompanying documentation (see ~/documentation/*). Notice that there must be\n\
#   only *one* source ON, otherwise the initialization procedure will issue an\n\
#   error message and stop the simulation.\n\
#\n\
[SECTION SOURCE BOX ISOTROPIC GAUSS SPECTRUM v.2006-08-01]\n\
(ON )                            STATUS (ON or OFF)\n\
 2                               PARTICLE TYPE (1=electron, 2=photon, 3=positron)\n\
 Energy(eV)  Probability         ENERGY SPECTRUM\n\
 25.0e3       1.0\n\
 25.0e3       0.000000010        # !! Force PENELOPE to create tables for high energy particles too\n\
  1.e6        -1\n";
	
	fprintf(peasy, "%s", text2);
	
	//center coordinates
	float cmCoorX = atof(xDim->value()) * 0.0001 / 2.0;
	float cmCoorY = atof(xDim->value()) * 0.0001 / 2.0;
	fprintf(peasy, " %.3e  %.3e  1.000e0     CENTER COORDINATES OF THE BOX ENCLOSURE (cm)\n", cmCoorX, cmCoorY);
	//fprintf("%.3e %.2e\n", cmCoorX, cmCoorY);
	 
	const char * text3 =
	" 3.000e-3  3.000e-3  0.000e0       BOX SIDES (cm)\n\
 0.000e0  0.000e0  0.000e0       EULER ANGLES [Rz,Ry,Rz](deg) TO ROTATE BOX\n\
 0                               MATERIAL (0=DO NOT CARE)\n\
 0.000e0  0.000e0  -1.000         DIRECTION VECTOR, NO NEED TO NORMALIZE\n\
 0.000e0                         ANGLE OF SEMIAPERTURE [0,180] (deg)\n\
[END OF BIGS SECTION]\n\
\n\
[SECTION SOURCE PHASE SPACE FILE v.2008-06-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0                               PSF FORMAT (0=STANDARD penEasy ASCII; 1=IAEA BINARY)\n\
(particles.psf                 ) PSF FILENAME (30 characters), REMOVE EXTENSION IF PSF FORMAT=1\n\
 1                               SPLITTING FACTOR\n\
 0.0  0.0  0.0                   EULER ANGLES [Rz,Ry,Rz](deg) TO ROTATE POSITION AND DIRECTION\n\
= 0.0  0.0  0.0                   CARTESIAN COMPONENTS [DX,DY,DZ](cm) OF POSITION SHIFT\n\
 1                               VALIDATE BEFORE SIMULATION (1=YES, MAY TAKE A WHILE; 0=NO)\n\
 0.000e0                         MAX PSF ENERGY (eV) (UNUSED IF VALIDATE=1; ADD 1023 keV FOR e+)\n\
[END OF SPSF SECTION]\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR SECTION PENGEOM+PENVOX\n\
#\n\
# * Three possible geometry models can be simulated: (i) quadrics; (ii) voxels\n\
#   (i.e. homogeneous parallelepipedic volume elements); and (iii) a mixture of\n\
#   quadrics and voxels.\n\
#\n\
# * For case (i), provide a file name in the QUADRICS field and leave the VOXELS\n\
#   field empty. The syntax of the quadrics file is described in the PENELOPE\n\
#   documentation.\n\
#\n\
# * For case (ii), do the reverse, leaving the QUADRICS field empty. The format\n\
#   of the voxels file is described in the file sample.vox included in this\n\
#   distribution. The voxels bounding box (i.e. the parallelepiped that delimits\n\
#   the set of voxels) is implicitly assumed to lie in the first octant of the\n\
#   simulation reference frame, that is, in the region {x>0,y>0,z>0}, with one\n\
#   of its corners located at the origin of coordinates. Its faces are therefore\n\
#   parallel to either the xy, xz or yz planes.\n\
#\n\
# * For case (iii), provide both file names. In this case quadric bodies are\n\
#   assumed to remove the voxels located \"underneath\", except for one body\n\
#   that we shall call the voxels body. The voxels body can be thought of as a\n\
#   \"hole\" in the quadric geometry through which voxels can be seen. The voxels\n\
#   body is identified by the user by giving the material of which it is made (in\n\
#   the field QUADRIC MATERIAL below). This material should not be present in any\n\
#   other body of the whole quadric geometry. Note also that the voxels body\n\
#   should not exceed in extension that of the parallelepipedic region where\n\
#   voxels are located, as defined in the voxels file.\n\
#\n\
# * For quadric+voxel geometries, each voxel mass is evaluated during the\n\
#   initialization by integrating density over voxel volume. The number of threads\n\
#   that are cast (along the z axis) through each voxel to perform this\n\
#   integration is determined by the value input in the GRANULARITY field. The\n\
#   actual number of threads equals the granularity squared. A reasonable value\n\
#   is 3; use larger values for spiky quadric geometries.\n\
#\n\
# * Note that vacuum bodies are NOT allowed inside the voxel region. This is to\n\
#   prevent possible miscalculations of a voxels mass.\n\
#\n\
[SECTION PENGEOM+PENVOX v.2008-06-01]\n\
(CsI.geo                       ) QUADRICS FILE NAME (30 characters), LEAVE EMPTY IF NONE\n\
(                              ) VOXELS FILE NAME (30 characters), LEAVE EMPTY IF NONE\n\
 1                               TRANSPARENT QUADRIC MAT (RELEVANT ONLY IF QUAD&VOX)\n\
 3                               GRANULARITY TO SCAN VOXELS (RELEVANT ONLY IF QUAD&VOX)\n\
[END OF GEO SECTION]\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR SECTION PENGEOM+PENVOX\n\
#\n\
# * The MATERIAL FILE must be a valid PENELOPE material file.\n\
#\n\
# * If the number of MATERIALS EXPECTED is set to zero, then this numbers is set\n\
#   automatically to the maximum material index found in the geometry files. The\n\
#   simulation parameters are then set, for all materials, as follows:\n\
#     - Eabs for electrons and positrons are both set to 1% of the initial source\n\
#       energy (E), with the limiting values of 50 eV (min) and 1 MeV (max).\n\
#     - Eabs for photons is set to 0.1% E with the limiting values of 50 eV and\n\
#       1 MeV.\n\
#     - C1 and C2 are both set to 0.1.\n\
#     - WCC is set to min(Eabs(e-),1% E)\n\
#     - WCR is set to min(Eabs(phot),0.1% E).\n\
#     - DSMAX is set to infinity.\n\
#   When this option is selected, the table header (MAT,EABS ... etc) and the\n\
#   succeeding lines with the values of these parameters must all be removed.\n\
#\n\
#\n\
#\n\
[SECTION PENELOPE v.2008-02-20]\n\
(CsI.mat                        ) MATERIAL FILE NAME (*** 30 characters ***)\n\
 1                               No. OF MATERIALS EXPECTED IN MAT FILE (0 for AUTO)\n\
 MAT  EABS(e-)  EABS(ph)  EABS(e+)  C1     C2     WCC       WCR       DSMAX   DESCRIPTION\n\
 1    50.0e0    50.0e0    50.00e0   0.01   0.01   50.00e0   50.00e0   0.020    CsI\n\
[END OF PEN SECTION]\n\
\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR THE TALLY SECTIONS\n\
#\n\
# * Details on the features and configuration of each tally are provided in their\n\
#   accompanying documentation (see ~/documentation/*.txt).\n\
#\n\
# * The required RELATIVE UNCERTAINTY that is specified for each tally (except for\n\
#   those that do not have an associated uncertainty, e.g. a phase-space file) is\n\
#   used as a condition to stop the simulation. Only when the requested\n\
#   relative uncertainties of *all* the tallies have been attained the uncertainty\n\
#   condition is considered fulfilled. Recall that the simulation can also be\n\
#   halted because the allotted time or the number of histories requested have\n\
#   been reached. Setting the RELATIVE UNCERTAINTY of all tallies to zero will\n\
#   prevent the execution from stopping for this cause.\n\
#\n\
# * Note for advanced users: when a certain tally scores nothing (i.e. zero) the\n\
#   corresponding REPORT routine reports 0% uncertainty but, at the same time, it\n\
#   reports that the requested uncertainty has not been reached, irrespective of\n\
#   the value introduced in the input file. This is to prevent the simulation from\n\
#   being stopped by a deceptive impression of accuracy in highly inefficient\n\
#   simulations, where the score and its standard deviation after a short period\n\
#   of time can be null.\n\
#\n\
[SECTION TALLY VOXEL DOSE v.2008-06-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0  0                            ROI MIN,MAX X-INDEX (0 0 FOR ALL VOXELS)\n\
 0  0                            ROI MIN,MAX Y-INDEX (0 0 FOR ALL VOXELS)\n\
 0  0                            ROI MIN,MAX Z-INDEX (0 0 FOR ALL VOXELS)\n\
 0                               INCLUDE QUAD. CONTRIBUTION TO VOXEL MASS & DOSE (1=YES,0=NO)\n\
 0                               PRINT VOXELS MASS IN REPORT (1=YES,0=NO)\n\
 0                               PRINT COORDINATES IN REPORT (1=YES,0=NO)\n\
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF VDD SECTION]\n\
\n\
[SECTION TALLY SPATIAL DOSE DISTRIB v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0.0  0.0     0                  XMIN,XMAX(cm),NXBIN (0 for DX=infty)\n\
 0.0  0.0     0                  YMIN,YMAX(cm),NYBIN (0 for DY=infty)\n\
-0.005 0.105  800                ZMIN,ZMAX(cm),NZBIN (0 for DZ=infty)\n\
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)\n\
 0.60                            RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF SDD SECTION]\n\
\n\
[SECTION TALLY CYLINDRICAL DOSE DISTRIB v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0.0 0.00001  500                RMIN,RMAX(cm),NRBIN (>0)\n\
 0.0 0.10     100                ZMIN,ZMAX(cm),NZBIN (0 for DZ=infty)\n\
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)\n\
 0.60                            RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF CDD SECTION]\n\
\n\
[SECTION TALLY SPHERICAL DOSE DISTRIB v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0.0 1.0   50                    RMIN,RMAX(cm),NRBIN (>0)\n\
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)\n\
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF SPD SECTION]\n\
\n\
[SECTION TALLY ENERGY DEPOSITION PULSE SPECTRUM v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 1                               DETECTION MATERIAL\n\
 0.0  1.0e9  100                 EMIN,EMAX(eV), No. OF E BINS\n\
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF EPS SECTION]\n\
\n\
[SECTION TALLY FLUENCE TRACK LENGTH v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 1                               DETECTION MATERIAL\n\
 1.0e2  1.0e9  70                EMIN,EMAX(eV), No. OF E BINS (LOG SCALE)\n\
 1.0e30                          RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF FTL SECTION]\n\
\n\
[SECTION TALLY PHASE SPACE FILE v.2008-06-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 0                               PSF FORMAT (0=STANDARD penEasy ASCII; 1=IAEA BINARY)\n\
 1                               DETECTION MATERIAL (NOT EQUAL 0)\n\
(output.psf                    ) PSF FILENAME (30 characters), REMOVE EXTENSION IF FORMAT=1\n\
[END OF PSF SECTION]\n\
\n\
[SECTION TALLY PARTICLE CURRENT SPECTRUM v.2006-08-01]\n\
(OFF)                            STATUS (ON or OFF)\n\
 1                               DETECTION MATERIAL\n\
 0.0 1.0e9   100                 EMIN,EMAX(eV), No. OF E BINS\n\
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED\n\
[END OF PCS SECTION]\n\
\n\
[SECTION TALLY PARTICLE TRACK STRUCTURE v.2008-05-15]\n\
(OFF)                            STATUS (ON or OFF)\n\
 100                             NUMBER OF HISTORIES TO DISPLAY (~100 RECOMMENDED)\n\
[END OF PTS SECTION]\n\
\n\
\n\
[SECTION TALLY ENERGY DEPOSITION EVENTS v.2010-06-03]\n\
(ON )                            STATUS (ON or OFF)\n\
 1                               DETECTOR SENSITIVE MATERIAL\n\
CsI_deposition_events.dat             NAME OUTPUT FILE\n\
[END OF EDE SECTION]\n\
\n\
\n\
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
# INSTRUCTIONS FOR THE INTERACTION FORCING SECTION\n\
#\n\
# * Interaction forcing is a variance reduction technique. It is described\n\
#   in detail in the PENELOPE manual. Use it judiciously.\n\
#\n\
# * Interaction forcing will only be applied if the particle's statistical\n\
#   weight (variable WGHT in common /TRACK/ of PENELOPE) falls in the\n\
#   specified interval [WMIN,WMAX], also known as the WEIGHT WINDOW.\n\
#   Otherwise, analog simulation is employed.\n\
#\n\
# * It is advisable to define a prudent weight window [WMIN,WMAX], such as\n\
#   [0.01,100], to prevent the occurrence of extreme statistical weights,\n\
#   which may give rise to low simulation efficiencies and numerical precision\n\
#   problems. Notice that a weight window such as [1.0,1.0] prevents the\n\
#   repetitive application of interaction forcing to secondary particles\n\
#   already generated in forced events.\n\
#\n\
# * One line must be entered for each combination of material and interaction\n\
#   type for which interaction forcing is to be applied. The contents of\n\
#   this line is:\n\
#   - Material number (MAT).\n\
#   - Particle type (KPAR). KPAR=1,2,3 is for electrons, photons and positrons,\n\
#     respectively.\n\
#   - Type of interaction (ICOL), see list of values below.\n\
#     By setting ICOL=0 all interactions are forced by the same amount.\n\
#   - Forcing factor (FORCING) by which the mean free path of interactions\n\
#     ICOL will be divided when particles of type KPAR are in material MAT.\n\
#\n\
# * PENELOPE labels the interaction mechanisms in the following way:\n\
#      Electrons (KPAR=1) and positrons (KPAR=3):\n\
#         ICOL = 1 artificial soft event (hinge).\n\
#              = 2 hard elastic collision.\n\
#              = 3 hard inelastic collision.\n\
#              = 4 hard bremsstrahlung emission.\n\
#              = 5 inner-shell ionization.\n\
#              = 6 positron annihilation.\n\
#              = 7 delta interaction.\n\
#              = 8 'auxiliary' fictitious interactions.\n\
#      Photons (KPAR=2):\n\
#         ICOL = 1 coherent (Rayleigh) scattering.\n\
#              = 2 incoherent (Compton) scattering.\n\
#              = 3 photoelectric absorption.\n\
#              = 4 electron-positron pair production.\n\
#              = 7 delta interaction.\n\
#              = 8 'auxiliary' fictitious interactions.\n\
#\n\
# * The last entered line (before END of SECTION) must have MAT=0, which\n\
#   signals the end of the list.\n\
#\n\
# * Beware of the fact that the use of interaction forcing may bias tallies\n\
#   based on pulse height spectra.\n\
#\n\
[SECTION INTERACTION FORCING v.2008-05-15]\n\
(OFF)                            STATUS (ON or OFF)\n\
 1.0  1.0                        WEIGHT WINDOW [WMIN,WMAX]\n\
 MAT  KPAR  ICOL  FORCING  (SET MAT=0 TO END LIST)\n\
 0    0     0     1.0\n\
[END OF VRIF SECTION]\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
# >>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";

	fprintf(peasy, "%s", text3);
	
	fclose(peasy);
}

/******************************
 * writeData
 * 
 * Saves updated inputs to
 * hybridMANTIS_input.in
 * ***************************/
void writeData()
{
	FILE * dataFP;
	
	dataFP = fopen("hybridMANTIS_input.in", "w");
	
	fprintf(dataFP,"# hybridMANTIS simulation input parameters\n");
	fprintf(dataFP,"#\n");
	fprintf(dataFP,"#	@file:		hybridMANTIS_input.in\n"); 
	fprintf(dataFP,"#	@author:	Diksha Sharma (Diksha.Sharma@fda.hhs.gov)\n");
	fprintf(dataFP,"#	@date:		Apr 9, 2012\n");
	fprintf(dataFP,"#\n\n");
	fprintf(dataFP,"%s     # number of x-ray histories to be simulated (N)\n", xRay->value());
	fprintf(dataFP,"%s          # min. number of optical photons that can be detected\n", minD->value());
	fprintf(dataFP,"%s       # max. (N*yield) number of optical photons that can be detected\n", maxD->value());
	fprintf(dataFP,"%s        # number of bins for storing pulse height specturm (maximum value=1000)\n", bins->value());
	fprintf(dataFP,"%s      # x-dimension of detector (in microns)\n", xDim->value());
	fprintf(dataFP,"%s      # y-dimension of detector (in microns)\n", yDim->value());
	fprintf(dataFP,"%s      # thickness of detector (in microns)\n", dThickness->value());
	fprintf(dataFP,"%s        # column radius (in microns)\n", colRad->value());
	fprintf(dataFP,"%s        # refractive index of column material\n", colRefrac->value());
	fprintf(dataFP,"%s        # refractive index of inter-columnar material\n", interColRefrac->value());
	fprintf(dataFP,"%s        # top surface absorption fraction\n", absorpFrac->value());
	fprintf(dataFP,"%s       # bulk absorption coefficient (in 1/microns)\n", bulkAbsorpCoeff->value());      
	fprintf(dataFP,"%s        # surface roughness coefficient\n", surRoughCoeff->value());       
	fprintf(dataFP,"%s        # minimum distance to the next column (in microns)\n", minDistance->value());
	fprintf(dataFP,"%s      # maximum distance to the next column (in microns)\n", maxDistance->value());
	fprintf(dataFP,"%s        # x lower bound of PRF image\n", xLower->value());       
	fprintf(dataFP,"%s        # y lower bound of PRF image\n", yLower->value());              
	fprintf(dataFP,"%s      # x upper bound of PRF image\n", xUpper->value());            
	fprintf(dataFP,"%s      # y upper bound of PRF image\n", yUpper->value());            
	fprintf(dataFP,"%s      # light yield (/eV)\n", lightYield->value());     
	fprintf(dataFP,"%s          # pixel pitch (in microns)  (max. pixels allowed in PRF image are 501x501. calculate this by {upper bound - lower bound}/pixel pitch.)\n", pixelPitch->value());       
	fprintf(dataFP,"%s       # non-ideal sensor reflectivity\n", niSensor->value());      
	fprintf(dataFP,"%s          # flag for running in the GPU (1) or only in the CPU (0)\n", flag->value());         
	fprintf(dataFP,"%s          # machine number\n", machine->value());
	 
	fclose(dataFP);
	
	dataFP = fopen("numPhotonHistories.txt", "w");
	fprintf(dataFP,"%s", numPhotonHistories->value());
	fclose(dataFP);
	
}
