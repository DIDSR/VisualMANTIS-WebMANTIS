<?php
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
//	Filename:	writeInputParameters.php
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Writes the user's input parameters into the correct files.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
$xRay = $_POST['xRay'];
$minDetect = $_POST['minDetect'];
$maxDetect = $_POST['maxDetect'];
$numBins = $_POST['numBins'];
$xDim = $_POST['xDim'];
$yDim = $_POST['yDim'];
$detThick = $_POST['detThick'];
$colRad = $_POST['colRad'];
$cri = $_POST['cri'];
$icri = $_POST['icri'];
$tsaf = $_POST['tsaf'];
$bac = $_POST['bac'];
$src = $_POST['src'];
$mindnc = $_POST['mindnc'];
$maxdnc = $_POST['maxdnc'];
$prfxlb = $_POST['prfxlb'];
$prfylb = $_POST['prfylb'];
$prfxub = $_POST['prfxub'];
$prfyub = $_POST['prfyub'];
$light = $_POST['light'];
$pixelPitch = $_POST['pixelPitch'];
$isr = $_POST['isr'];
$histories = $_POST['histories'];
$folder = $_POST['folder'];

//write to numPhotonHistories
$fp2 = fopen($folder."/numPhotonHistories.txt", 'w');
fwrite($fp2, $histories."\n");
fclose($fp2);

//write to hybridMANTIS_input.in
$fp = fopen($folder."/hybridMANTIS_input.in", 'w');
fwrite($fp, "# hybridMANTIS simulation input parameters\n");
fwrite($fp, "#\n");
fwrite($fp, "#	@file:		hybridMANTIS_input.in\n");
fwrite($fp, "#	@author:	Diksha Sharma (Diksha.Sharma@fda.hhs.gov)\n");
fwrite($fp, "#	@date:		Apr 9, 2012\n");
fwrite($fp, "#\n");
fwrite($fp, "\n");

fwrite($fp, $xRay."\t# number of x-ray histories to be simulated (N)\n");
fwrite($fp, $minDetect."\t# min. number of optical photons that can be detected\n");
fwrite($fp, $maxDetect."\t# max. (N*yield) number of optical photons that can be detected\n");
fwrite($fp, $numBins."\t# number of bins for storing pulse height specturm (maximum value=1000)\n");
fwrite($fp, $xDim."\t# x-dimension of detector (in microns)\n");
fwrite($fp, $yDim."\t# y-dimension of detector (in microns)\n");
fwrite($fp, $detThick."\t# thickness of detector (in microns)\n");
fwrite($fp, $colRad."\t# column radius (in microns)\n");
fwrite($fp, $cri."\t# refractive index of column material\n");
fwrite($fp, $icri."\t# refractive index of inter-columnar material\n");
fwrite($fp, $tsaf."\t# top surface absorption fraction\n");
fwrite($fp, $bac."\t# bulk absorption coefficient (in 1/microns)\n");
fwrite($fp, $src."\t# surface roughness coefficient\n");
fwrite($fp, $mindnc."\t# minimum distance to the next column (in microns)\n");
fwrite($fp, $maxdnc."\t# maximum distance to the next column (in microns)\n");
fwrite($fp, $prfxlb."\t# x lower bound of PRF image\n");
fwrite($fp, $prfylb."\t# y lower bound of PRF image\n");
fwrite($fp, $prfxub."\t# x upper bound of PRF image\n");
fwrite($fp, $prfyub."\t# y upper bound of PRF image\n");
fwrite($fp, $light."\t# light yield (/eV)\n");
fwrite($fp, $pixelPitch."\t# pixel pitch (in microns)  (max. pixels allowed in PRF image are 501x501. calculate this by {upper bound - lower bound}/pixel pitch.)\n");
fwrite($fp, $isr."\t# non-ideal sensor reflectivity\n");
fwrite($fp, "1\t# flag for running in the GPU (1) or only in the CPU (0)\n");
fwrite($fp, "1\t# machine number\n");
//fwrite($fp, $histories."\n");
fclose($fp);

//write to saved.txt
$fp = fopen($folder."/saved.txt", 'w');
fwrite($fp, $xRay."\n");
fwrite($fp, $minDetect."\n");
fwrite($fp, $maxDetect."\n");
fwrite($fp, $numBins."\n");
fwrite($fp, $xDim."\n");
fwrite($fp, $yDim."\n");
fwrite($fp, $detThick."\n");
fwrite($fp, $colRad."\n");
fwrite($fp, $cri."\n");
fwrite($fp, $icri."\n");
fwrite($fp, $tsaf."\n");
fwrite($fp, $bac."\n");
fwrite($fp, $src."\n");
fwrite($fp, $mindnc."\n");
fwrite($fp, $maxdnc."\n");
fwrite($fp, $prfxlb."\n");
fwrite($fp, $prfylb."\n");
fwrite($fp, $prfxub."\n");
fwrite($fp, $prfyub."\n");
fwrite($fp, $light."\n");
fwrite($fp, $pixelPitch."\n");
fwrite($fp, $isr."\n");
fwrite($fp, "1\n");
fwrite($fp, "1\n");
fwrite($fp, $histories."\n");
fclose($fp);

//write to CsI.geo
/*$fp = fopen($folder."/CsI.geo", 'w');

fwrite($fp, 
">>>> PenEasy geometry: CsI rectangle in vacuum >>>>



                      ^
                     /|\ z-axis
                      |
             vacuum   |
                      |
                      |               
      111111111111111111111111111111111 -> x=T1 cm  Cesium iodide
      111111111111111111111111111111111
      111111111111111111111111111111111---->  x=0 cm

The actual geometry definition begins here:
0000000000000000000000000000000000000000000000000000000000000000\n");

$det = (int)$detThick / 2; 

fwrite($fp, "SURFACE (   1)   Plane Z = ".(-1 * $det)." um (".(-1 * (float)($det / 10000))." cm)
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-7.500000000000000E-03,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)   Plane Z = ".($det)." um (".(float)($det / 10000)." cm)
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(+7.500000000000000E-03,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)   Plane X = 0 cm
INDICES=( 0, 0, 0, 0, 0)
     AX=(+1.000000000000000E+00,   0)
     A0=(+0.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)   Plane X = 909 um (0.0909 cm)
INDICES=( 0, 0, 0, 0, 0)
     AX=(+1.000000000000000E+00,   0)
     A0=(-9.090000000000000E-02,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   5)   Plane Y = 0 cm
INDICES=( 0, 0, 0, 0, 0)
     AY=(+1.000000000000000E+00,   0)
     A0=(+0.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   6)   Plane Y = 909 um (0.0909 cm)
INDICES=( 0, 0, 0, 0, 0)
     AY=(+1.000000000000000E+00,   0)
     A0=(-9.090000000000000E-02,   0)
0000000000000000000000000000000000000000000000000000000000000000\n");

fwrite($fp, "
BODY    (   1)   CsI detector
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=(+1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=(+1)
SURFACE (   4), SIDE POINTER=(-1)
SURFACE (   5), SIDE POINTER=(+1)
SURFACE (   6), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000

>>> END OF INPUT DATA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





GEOMETRY FILE LAYOUT:

0000000000000000000000000000000000000000000000000000000000000000
SURFACE (    )   reduced form
INDICES=( 1, 1, 1, 1, 1)
X-SCALE=(+1.000000000000000E+00,   0)              (DEFAULT=1.0)
Y-SCALE=(+1.000000000000000E+00,   0)              (DEFAULT=1.0)
Z-SCALE=(+1.000000000000000E+00,   0)              (DEFAULT=1.0)
  OMEGA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
  THETA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
    PHI=(+0.000000000000000E+00,   0) RAD          (DEFAULT=0.0)
X-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Y-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Z-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (    )   implicit form
INDICES=( 0, 0, 0, 0, 0)
    AXX=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
    AXY=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
    AXZ=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
    AYY=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
    AYZ=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
    AZZ=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
     AX=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
     AY=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
     AZ=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
     A0=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
1111111111111111111111111111111111111111111111111111111111111111
  OMEGA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
  THETA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
    PHI=(+0.000000000000000E+00,   0) RAD          (DEFAULT=0.0)
X-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Y-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Z-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (    )   text
MATERIAL(    )
SURFACE (    ), SIDE POINTER=( 1)
BODY    (    )
MODULE  (    )
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (    )   text
MATERIAL(    )
SURFACE (    ), SIDE POINTER=( 1)
BODY    (    )
MODULE  (    )
1111111111111111111111111111111111111111111111111111111111111111
  OMEGA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
  THETA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
    PHI=(+0.000000000000000E+00,   0) RAD          (DEFAULT=0.0)
X-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Y-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Z-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
CLONE   (    )   copies one module and moves it
MODULE  (    )   original module
1111111111111111111111111111111111111111111111111111111111111111
  OMEGA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
  THETA=(+0.000000000000000E+00,   0) DEG          (DEFAULT=0.0)
    PHI=(+0.000000000000000E+00,   0) RAD          (DEFAULT=0.0)
X-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Y-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
Z-SHIFT=(+0.000000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
INCLUDE
   FILE=(filename.ext)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000


>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

");

fclose($fp);*/

//write to penEasy_CsI_input.in
$fp = fopen($folder."/penEasy_CsI_input.in", 'w');
fwrite($fp,	"# >>>> INPUT FILE FOR penEasy >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# CASE DESCRIPTION:
#   penEasy with No EM fields!
#
# LAST UPDATE:
#   2008-04-19 by JS
#   2011-06-06 by DS


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR SECTION CONFIG
#
# * The simulation will be finished when ANY of the following conditions is
#   fulfilled: i) the requested number of histories has been reached; ii) the
#   alloted time has been exhausted; or iii) ALL the requested relative
#   uncertainties have been reached.
#
# * Allotted time is interpreted as real time if it is a positive number.
#   Otherwise, CPU (i.e. user) time is assumed.
#
# * Similarly, the update period is interpreted as real time if it is positive;
#   otherwise it is the number of histories between updates. See the README file
#   for a description of the tasks performed during each update. Note that the
#   real time between updates cannot be larger than 50000 s.
#
# * Setting both random seeds equal to zero is an indication that they must be read
#   from an external file. The name of the file containing the seeds is read from
#   an additional input line that must be introduced after that of 'RANDOM SEEDS'.
#   This new line must have the following format:
#(rngseed.in                    ) RANDOM SEEDS FILE NAME (*** 30 characters ***)
#   In this example, the file 'rngseed.in' is opened and read. The two seeds must
#   be on the first line of rngseed.in, separated by one or more blanks. This
#   feature is useful when using parallel computing (see README file).
#
[SECTION CONFIG v.2006-08-01]\n");

$svar = (string)$xRay;
fwrite($fp, "  ".$svar[0].".".substr($svar, 1)."e".strlen(substr($svar, 1))."                            NO. OF HISTORIES (<1e15)\n");

fwrite($fp, "  5.0e15                           ALLOTTED TIME (s) (+ FOR REAL TIME; - FOR CPU TIME)
  50000.0                          UPDATE INTERVAL (+ FOR REAL TIME (s) < 50000; - FOR HISTORIES)
 3141592   20090625                INITIAL RANDOM SEEDS
[END OF CONFIG SECTION]


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR SOURCE SECTIONS
#
# * Details on the features and configuration of each source are provided in their
#   accompanying documentation (see ~/documentation/*). Notice that there must be
#   only *one* source ON, otherwise the initialization procedure will issue an
#   error message and stop the simulation.
#
[SECTION SOURCE BOX ISOTROPIC GAUSS SPECTRUM v.2006-08-01]
(ON )                            STATUS (ON or OFF)
 2                               PARTICLE TYPE (1=electron, 2=photon, 3=positron)
 Energy(eV)  Probability         ENERGY SPECTRUM
 25.0e3       1.0
 25.0e3       0.000000010        # !! Force PENELOPE to create tables for high energy particles too
  1.e6        -1\n");

//center coordinates
$cmCoorX = (float)$xDim * 0.0001 / 2.0;
$cmCoorY = (float)$yDim * 0.0001 / 2.0;

fprintf($fp, " %.03e  %.03e  1.000e0     CENTER COORDINATES OF THE BOX ENCLOSURE (cm)\n", $cmCoorX, $cmCoorY);
fwrite($fp, " 3.000e-3  3.000e-3  0.000e0       BOX SIDES (cm)
 0.000e0  0.000e0  0.000e0       EULER ANGLES [Rz,Ry,Rz](deg) TO ROTATE BOX
 0                               MATERIAL (0=DO NOT CARE)
 0.000e0  0.000e0  -1.000         DIRECTION VECTOR, NO NEED TO NORMALIZE
 0.000e0                         ANGLE OF SEMIAPERTURE [0,180] (deg)
[END OF BIGS SECTION]

[SECTION SOURCE PHASE SPACE FILE v.2008-06-01]
(OFF)                            STATUS (ON or OFF)
 0                               PSF FORMAT (0=STANDARD penEasy ASCII; 1=IAEA BINARY)
(particles.psf                 ) PSF FILENAME (30 characters), REMOVE EXTENSION IF PSF FORMAT=1
 1                               SPLITTING FACTOR
 0.0  0.0  0.0                   EULER ANGLES [Rz,Ry,Rz](deg) TO ROTATE POSITION AND DIRECTION
= 0.0  0.0  0.0                   CARTESIAN COMPONENTS [DX,DY,DZ](cm) OF POSITION SHIFT
 1                               VALIDATE BEFORE SIMULATION (1=YES, MAY TAKE A WHILE; 0=NO)
 0.000e0                         MAX PSF ENERGY (eV) (UNUSED IF VALIDATE=1; ADD 1023 keV FOR e+)
[END OF SPSF SECTION]


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR SECTION PENGEOM+PENVOX
#
# * Three possible geometry models can be simulated: (i) quadrics; (ii) voxels
#   (i.e. homogeneous parallelepipedic volume elements); and (iii) a mixture of
#   quadrics and voxels.
#
# * For case (i), provide a file name in the QUADRICS field and leave the VOXELS
#   field empty. The syntax of the quadrics file is described in the PENELOPE
#   documentation.
#
# * For case (ii), do the reverse, leaving the QUADRICS field empty. The format
#   of the voxels file is described in the file sample.vox included in this
#   distribution. The voxels bounding box (i.e. the parallelepiped that delimits
#   the set of voxels) is implicitly assumed to lie in the first octant of the
#   simulation reference frame, that is, in the region {x>0,y>0,z>0}, with one
#   of its corners located at the origin of coordinates. Its faces are therefore
#   parallel to either the xy, xz or yz planes.
#
# * For case (iii), provide both file names. In this case quadric bodies are
#   assumed to remove the voxels located \"underneath\", except for one body
#   that we shall call the voxels body. The voxels body can be thought of as a
#   \"hole\" in the quadric geometry through which voxels can be seen. The voxels
#   body is identified by the user by giving the material of which it is made (in
#   the field QUADRIC MATERIAL below). This material should not be present in any
#   other body of the whole quadric geometry. Note also that the voxels body
#   should not exceed in extension that of the parallelepipedic region where
#   voxels are located, as defined in the voxels file.
#
# * For quadric+voxel geometries, each voxel mass is evaluated during the
#   initialization by integrating density over voxel volume. The number of threads
#   that are cast (along the z axis) through each voxel to perform this
#   integration is determined by the value input in the GRANULARITY field. The
#   actual number of threads equals the granularity squared. A reasonable value
#   is 3; use larger values for spiky quadric geometries.
#
# * Note that vacuum bodies are NOT allowed inside the voxel region. This is to
#   prevent possible miscalculations of a voxels mass.
#
[SECTION PENGEOM+PENVOX v.2008-06-01]
(CsI.geo                       ) QUADRICS FILE NAME (30 characters), LEAVE EMPTY IF NONE
(                              ) VOXELS FILE NAME (30 characters), LEAVE EMPTY IF NONE
 1                               TRANSPARENT QUADRIC MAT (RELEVANT ONLY IF QUAD&VOX)
 3                               GRANULARITY TO SCAN VOXELS (RELEVANT ONLY IF QUAD&VOX)
[END OF GEO SECTION]


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR SECTION PENGEOM+PENVOX
#
# * The MATERIAL FILE must be a valid PENELOPE material file.
#
# * If the number of MATERIALS EXPECTED is set to zero, then this numbers is set
#   automatically to the maximum material index found in the geometry files. The
#   simulation parameters are then set, for all materials, as follows:
#     - Eabs for electrons and positrons are both set to 1% of the initial source
#       energy (E), with the limiting values of 50 eV (min) and 1 MeV (max).
#     - Eabs for photons is set to 0.1% E with the limiting values of 50 eV and
#       1 MeV.
#     - C1 and C2 are both set to 0.1.
#     - WCC is set to min(Eabs(e-),1% E)
#     - WCR is set to min(Eabs(phot),0.1% E).
#     - DSMAX is set to infinity.
#   When this option is selected, the table header (MAT,EABS ... etc) and the
#   succeeding lines with the values of these parameters must all be removed.
#
#
#
[SECTION PENELOPE v.2008-02-20]
(CsI.mat                        ) MATERIAL FILE NAME (*** 30 characters ***)
 1                               No. OF MATERIALS EXPECTED IN MAT FILE (0 for AUTO)
 MAT  EABS(e-)  EABS(ph)  EABS(e+)  C1     C2     WCC       WCR       DSMAX   DESCRIPTION
 1    50.0e0    50.0e0    50.00e0   0.01   0.01   50.00e0   50.00e0   0.020    CsI
[END OF PEN SECTION]



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR THE TALLY SECTIONS
#
# * Details on the features and configuration of each tally are provided in their
#   accompanying documentation (see ~/documentation/*.txt).
#
# * The required RELATIVE UNCERTAINTY that is specified for each tally (except for
#   those that do not have an associated uncertainty, e.g. a phase-space file) is
#   used as a condition to stop the simulation. Only when the requested
#   relative uncertainties of *all* the tallies have been attained the uncertainty
#   condition is considered fulfilled. Recall that the simulation can also be
#   halted because the allotted time or the number of histories requested have
#   been reached. Setting the RELATIVE UNCERTAINTY of all tallies to zero will
#   prevent the execution from stopping for this cause.
#
# * Note for advanced users: when a certain tally scores nothing (i.e. zero) the
#   corresponding REPORT routine reports 0% uncertainty but, at the same time, it
#   reports that the requested uncertainty has not been reached, irrespective of
#   the value introduced in the input file. This is to prevent the simulation from
#   being stopped by a deceptive impression of accuracy in highly inefficient
#   simulations, where the score and its standard deviation after a short period
#   of time can be null.
#
[SECTION TALLY VOXEL DOSE v.2008-06-01]
(OFF)                            STATUS (ON or OFF)
 0  0                            ROI MIN,MAX X-INDEX (0 0 FOR ALL VOXELS)
 0  0                            ROI MIN,MAX Y-INDEX (0 0 FOR ALL VOXELS)
 0  0                            ROI MIN,MAX Z-INDEX (0 0 FOR ALL VOXELS)
 0                               INCLUDE QUAD. CONTRIBUTION TO VOXEL MASS & DOSE (1=YES,0=NO)
 0                               PRINT VOXELS MASS IN REPORT (1=YES,0=NO)
 0                               PRINT COORDINATES IN REPORT (1=YES,0=NO)
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED
[END OF VDD SECTION]

[SECTION TALLY SPATIAL DOSE DISTRIB v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 0.0  0.0     0                  XMIN,XMAX(cm),NXBIN (0 for DX=infty)
 0.0  0.0     0                  YMIN,YMAX(cm),NYBIN (0 for DY=infty)
-0.005 0.105  800                ZMIN,ZMAX(cm),NZBIN (0 for DZ=infty)
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)
 0.60                            RELATIVE UNCERTAINTY (%) REQUESTED
[END OF SDD SECTION]

[SECTION TALLY CYLINDRICAL DOSE DISTRIB v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 0.0 0.00001  500                RMIN,RMAX(cm),NRBIN (>0)
 0.0 0.10     100                ZMIN,ZMAX(cm),NZBIN (0 for DZ=infty)
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)
 0.60                            RELATIVE UNCERTAINTY (%) REQUESTED
[END OF CDD SECTION]

[SECTION TALLY SPHERICAL DOSE DISTRIB v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 0.0 1.0   50                    RMIN,RMAX(cm),NRBIN (>0)
 1                               PRINT COORDINATES IN REPORT (1=YES,0=NO)
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED
[END OF SPD SECTION]

[SECTION TALLY ENERGY DEPOSITION PULSE SPECTRUM v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 1                               DETECTION MATERIAL
 0.0  1.0e9  100                 EMIN,EMAX(eV), No. OF E BINS
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED
[END OF EPS SECTION]

[SECTION TALLY FLUENCE TRACK LENGTH v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 1                               DETECTION MATERIAL
 1.0e2  1.0e9  70                EMIN,EMAX(eV), No. OF E BINS (LOG SCALE)
 1.0e30                          RELATIVE UNCERTAINTY (%) REQUESTED
[END OF FTL SECTION]

[SECTION TALLY PHASE SPACE FILE v.2008-06-01]
(OFF)                            STATUS (ON or OFF)
 0                               PSF FORMAT (0=STANDARD penEasy ASCII; 1=IAEA BINARY)
 1                               DETECTION MATERIAL (NOT EQUAL 0)
(output.psf                    ) PSF FILENAME (30 characters), REMOVE EXTENSION IF FORMAT=1
[END OF PSF SECTION]

[SECTION TALLY PARTICLE CURRENT SPECTRUM v.2006-08-01]
(OFF)                            STATUS (ON or OFF)
 1                               DETECTION MATERIAL
 0.0 1.0e9   100                 EMIN,EMAX(eV), No. OF E BINS
 0.0                             RELATIVE UNCERTAINTY (%) REQUESTED
[END OF PCS SECTION]

[SECTION TALLY PARTICLE TRACK STRUCTURE v.2008-05-15]
(OFF)                            STATUS (ON or OFF)
 100                             NUMBER OF HISTORIES TO DISPLAY (~100 RECOMMENDED)
[END OF PTS SECTION]


[SECTION TALLY ENERGY DEPOSITION EVENTS v.2010-06-03]
(ON )                            STATUS (ON or OFF)
 1                               DETECTOR SENSITIVE MATERIAL
CsI_deposition_events.dat             NAME OUTPUT FILE
[END OF EDE SECTION]


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INSTRUCTIONS FOR THE INTERACTION FORCING SECTION
#
# * Interaction forcing is a variance reduction technique. It is described
#   in detail in the PENELOPE manual. Use it judiciously.
#
# * Interaction forcing will only be applied if the particle's statistical
#   weight (variable WGHT in common /TRACK/ of PENELOPE) falls in the
#   specified interval [WMIN,WMAX], also known as the WEIGHT WINDOW.
#   Otherwise, analog simulation is employed.
#
# * It is advisable to define a prudent weight window [WMIN,WMAX], such as
#   [0.01,100], to prevent the occurrence of extreme statistical weights,
#   which may give rise to low simulation efficiencies and numerical precision
#   problems. Notice that a weight window such as [1.0,1.0] prevents the
#   repetitive application of interaction forcing to secondary particles
#   already generated in forced events.
#
# * One line must be entered for each combination of material and interaction
#   type for which interaction forcing is to be applied. The contents of
#   this line is:
#   - Material number (MAT).
#   - Particle type (KPAR). KPAR=1,2,3 is for electrons, photons and positrons,
#     respectively.
#   - Type of interaction (ICOL), see list of values below.
#     By setting ICOL=0 all interactions are forced by the same amount.
#   - Forcing factor (FORCING) by which the mean free path of interactions
#     ICOL will be divided when particles of type KPAR are in material MAT.
#
# * PENELOPE labels the interaction mechanisms in the following way:
#      Electrons (KPAR=1) and positrons (KPAR=3):
#         ICOL = 1 artificial soft event (hinge).
#              = 2 hard elastic collision.
#              = 3 hard inelastic collision.
#              = 4 hard bremsstrahlung emission.
#              = 5 inner-shell ionization.
#              = 6 positron annihilation.
#              = 7 delta interaction.
#              = 8 'auxiliary' fictitious interactions.
#      Photons (KPAR=2):
#         ICOL = 1 coherent (Rayleigh) scattering.
#              = 2 incoherent (Compton) scattering.
#              = 3 photoelectric absorption.
#              = 4 electron-positron pair production.
#              = 7 delta interaction.
#              = 8 'auxiliary' fictitious interactions.
#
# * The last entered line (before END of SECTION) must have MAT=0, which
#   signals the end of the list.
#
# * Beware of the fact that the use of interaction forcing may bias tallies
#   based on pulse height spectra.
#
[SECTION INTERACTION FORCING v.2008-05-15]
(OFF)                            STATUS (ON or OFF)
 1.0  1.0                        WEIGHT WINDOW [WMIN,WMAX]
 MAT  KPAR  ICOL  FORCING  (SET MAT=0 TO END LIST)
 0    0     0     1.0
[END OF VRIF SECTION]





























# >>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");	
fclose($fp);

$ret = exec("cp ".$folder."/*.in ".$folder."/visualmantis/");
$ret = exec("cp ".$folder."/*.txt ".$folder."/visualmantis/");
$ret = exec("cp ".$folder."/*.geo ".$folder."/visualmantis/");

$response = array();
$response['ret'] = 'true';
echo json_encode($response);
 
?>
