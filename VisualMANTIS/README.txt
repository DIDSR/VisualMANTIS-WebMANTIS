##########################################################################################################################################################
#
# ****Disclaimer****
#  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
#  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
#  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
#  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
#  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
#  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
#  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
#  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, 
# and any modified versions bear some notice that they have been modified.
#
#	@file    README.txt
#   @author  Han Dong, Diksha Sharma (Diksha.Sharma@fda.hhs.gov), Aldo Badano (Aldo.Badano@fda.hhs.gov)
#   @date    Apr 19, 2012
#
##########################################################################################################################################################


*************
INTRODUCTION
*************

visualMANTIS is a visualization tool that is built on top of hybridMANTIS - a Monte Carlo tool for modeling x-ray detectors with columnar scintillators.

The source code is free and open software in the public domain, as explained in the Disclaimer section above. 
The software distribution website is: https://github.com/diamfda


*******************************
CODE COMPILATION AND EXECUTION
*******************************

visualMANTIS has been tested only on the Linux operating system. The CUDA libraries, GNU gcc and gfortran compiler and GNU scientific library needs to be pre installed before running. visualMANTIS relies on having a prior version of hybridmantis and needs to be combined with the files in hybridMANTIS . A bash script is included to compile the CUDA and C codes with PENELOPE and penEasy files. This file may have to be edited to modify the library paths. A pre-compiled executable is also attached. It was compiled using CUDA version 4.2, gcc version 4.4.5, gsl version 1.14, fltk version 1.3.1, and gnuplot 4.6.

To run visualMANTIS:
	./visualmantisGUI
	
Compilation steps are listed in detail in 'MANUAL_visualMANTIS.pdf'.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																//
//    NOTE: PENELOPE 2006 SOURCE CODE FILES PENELOPE.F AND PENGEOM.F ARE NOT DISTRIBUTED WITH hybridMANTIS PACKAGE, 		//
//    BUT ARE NEEDED FOR COMPILING. IF THESE FILES ARE NEEDED, PLEASE CONTACT DIKSHA SHARMA AT Diksha.Sharma(at)fda.hhs.gov.	//
//																//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


**********************************
hybridMANTIS v1.0 PACKAGE CONTENTS
**********************************

	* modified visualMANTIS kernel and C code:  'visualMANTIS_****' CUDA and C files. 
	
	* visualization : contains '****.cxx' files, these are the FLTK and C code to handle the GUI.
	
	* visualizationGUI:	 		main GUI file to launch visualMANTIS
	
	* visualMANTIS_ver1_0.x: 	visualMANTIS v1.0 executable - this only handles retrieving visualization data from CUDA kernels
	
	* compile_visualmantis_ver1_0.sh:		Script for compiling visualMANTIS. The user will require two more files PENELOPE.F, PENGEOM.F to compile this package successfully. These files are not distributed with this package. If these files are needed, please contact Diksha Sharma at Diksha.Sharma (at) fda.hhs.gov.
	
	* example/ : contains sample input and output files
	
	* README.txt:			this file.
	

For more details, read 'MANUAL_visualMANTIS.pdf'.
	
[1] Sharma, D.; Badal, A.; and Badano, A. 2012. hybridMANTIS: a CPU–GPU Monte Carlo method for modeling indirect x-ray detectors with columnar scintillators. Physics in Medicine and Biology 57(8):2357.
