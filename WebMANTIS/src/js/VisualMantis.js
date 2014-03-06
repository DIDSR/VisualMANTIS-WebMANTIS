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
//	Filename:	VisualMantis.js
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Contains global variables needed
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
var currentUser = '';
var fetchCounter = 0;
var phsCounter = 1;
var prfCounter = 1;
var currentXray = 0;
var globalJSONData = {};
window.jsonData = {};
var myPID = 0;
var openvar = 0;
var objects = new Array();
var spheres = new Array();
var ceilings = new Array();
var floors = new Array();
var lines = new Array();
var prfxub2, prfyub2, photons2, xray2, mindnc2, maxdnc2, numBins2, xdim2, ydim2, detector2, colRad2, cri2, icri2, tsaf2, bac2, src2, prfxlb2, prfylb2, light2, pixel2, nisr2, mind2, maxd2 = 0;
