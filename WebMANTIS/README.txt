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
#   @author  Han Dong, Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#   @date    2/27/2014
#
##########################################################################################################################################################


*************
INTRODUCTION
*************

WebMANTIS is a graphical interface that uses VisualMANTIS and hybridMANTIS (a Monte Carlo tool for modeling x-ray detectors with columnar scintillators) [1] to display simulation information and visualizations in a web browser. WebMANTIS borrows heavily from utility files and sample codes that were provided by WebGL, more information about WebGL can be found here (http://threejs.org/).

The source code is free and open software in the public domain, as explained in the Disclaimer section above. 
The software distribution website is: https://github.com/diamfda. 


*******************************
CODE COMPILATION AND EXECUTION
*******************************

WebMANTIS has been tested only on the Ubuntu operating system, it requires the installation of the Ubuntu LAMP (ApacheMySQLPHP) package in order to set up the local webserver. The code doesn't need to be compiled since it is built using scripting languages, after Ubuntu LAMP has been set up this code can then copied into the webserver directory and it will be accessible via the web.
	
Compilation steps and libraries required are listed in detail in 'MANUAL_WebMANTIS.pdf'.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																//
//    NOTE: PENELOPE 2006 SOURCE CODE FILES PENELOPE.F AND PENGEOM.F ARE NOT DISTRIBUTED WITH hybridMANTIS PACKAGE, 		//
//    BUT ARE NEEDED FOR COMPILING. IF THESE FILES ARE NEEDED, PLEASE CONTACT DIKSHA SHARMA AT Diksha.Sharma(at)fda.hhs.gov.	//
//																//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


**********************************
WebMANTIS v1.0 PACKAGE CONTENTS
**********************************

	* build/ : contains the WebGL JavaScript files needed for rendering the 3D objects
	
	* examples/ : contains utility libraries that came with WebGL and WebMANTIS also utilizes them
	
	* src/ : the main PHP and Python backend scripts are located here to communicate the results of the GPU running hybridMANTIS to the front end WebGL visualization
	
	* src/Guest : an example of a sample user using WebMANTIS, the contents of Guest contains necessary simulation input files and the VisualMANTIS component that utilizes the GPU to run the simulation
	
	* src/js/ : JavaScript files that form the main graphical interface, controls different parts of the interface such as the  menubar, sidebar, toolbar and main 3D view
	
	* src/libs/ : utility JavaScript libraries provided by WebGL package
	
	* sample_output/ : contains examples of all the different types of output produced by WebMANTIS
	
	* libevent-multithreaded-server/ : contains third party libevent code that acts the intermediate server
	
	* README.txt : provides information on  various parts of the package contents
	

For more details, read 'MANUAL_WebMANTIS.pdf'.
	
[1] Sharma, D.; Badal, A.; and Badano, A. 2012. hybridMANTIS: a CPU–GPU Monte Carlo method for modeling indirect x-ray detectors with columnar scintillators. Physics in Medicine and Biology 57(8):2357.
