
INTRODUCTION
=============

WebMANTIS is a graphical interface that uses VisualMANTIS and hybridMANTIS (a Monte Carlo tool for modeling x-ray detectors with columnar scintillators) to display simulation information and visualizations in a web browser. WebMANTIS borrows heavily from utility files and sample codes that were provided by WebGL, more information about WebGL can be found here (http://threejs.org/).

CODE COMPILATION AND EXECUTION
===============================

WebMANTIS has been tested only on the Ubuntu operating system, it requires the installation of the Ubuntu LAMP (ApacheMySQLPHP) package in order to set up the local webserver. The code doesn't need to be compiled since it is built using scripting languages, after Ubuntu LAMP has been set up this code can then copied into the webserver directory and it will be accessible via the web.
	
Compilation steps and libraries required are listed in detail in 'VisualMANTIS_WebMANTIS_manual.pdf'.

NOTE: PENELOPE 2006 SOURCE CODE FILES PENELOPE.F AND PENGEOM.F ARE NOT DISTRIBUTED WITH hybridMANTIS PACKAGE, BUT ARE NEEDED FOR COMPILING. IF THESE FILES ARE NEEDED, PLEASE CONTACT DIKSHA SHARMA AT Diksha.Sharma(at)fda.hhs.gov.

WebMANTIS v1.0 PACKAGE CONTENTS
================================

	* build/ : contains the WebGL JavaScript files needed for rendering the 3D objects
	
	* examples/ : contains utility libraries that came with WebGL and WebMANTIS also utilizes them
	
	* src/ : the main PHP and Python backend scripts are located here to communicate the results of the GPU running hybridMANTIS to the front end WebGL visualization
	
	* src/Guest : an example of a sample user using WebMANTIS, the contents of Guest contains necessary simulation input files and the VisualMANTIS component that utilizes the GPU to run the simulation
	
	* src/js/ : JavaScript files that form the main graphical interface, controls different parts of the interface such as the  menubar, sidebar, toolbar and main 3D view
	
	* src/libs/ : utility JavaScript libraries provided by WebGL package
	
	* sample_output/ : contains examples of all the different types of output produced by WebMANTIS
	
	* libevent-multithreaded-server/ : contains third party libevent code that acts the intermediate server
	
	* README.md : provides information on  various parts of the package contents
