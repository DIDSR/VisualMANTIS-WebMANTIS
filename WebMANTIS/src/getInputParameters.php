<?
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
//	Filename:	getInputParameters.php
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This code saves the input arguments for future use 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
 
$folder = $_POST['folder'];
$f = fopen($folder."/saved.txt", "r"); 

$response = array();
	
if($f)
{
	$response['ret'] = "True";
	$response['xray'] = trim(fgets($f, 4096));
	$response['minD'] = trim(fgets($f, 4096));
	$response['maxD'] = trim(fgets($f, 4096));
	$response['numBins'] = trim(fgets($f, 4096));
	$response['xDim'] = trim(fgets($f, 4096));
	$response['yDim'] = trim(fgets($f, 4096));
	$response['detector'] = trim(fgets($f, 4096));
	$response['colRad'] = trim(fgets($f, 4096));
	$response['cri'] = trim(fgets($f, 4096));
	$response['icri'] = trim(fgets($f, 4096));
	$response['tsaf'] = trim(fgets($f, 4096));
	$response['bac'] = trim(fgets($f, 4096));
	$response['src'] = trim(fgets($f, 4096));
	$response['mindnc'] = trim(fgets($f, 4096));
	$response['maxdnc'] = trim(fgets($f, 4096));
	$response['prfxlb'] = trim(fgets($f, 4096));
	$response['prfylb'] = trim(fgets($f, 4096));
	$response['prfxub'] = trim(fgets($f, 4096));
	$response['prfyub'] = trim(fgets($f, 4096));
	$response['light'] = trim(fgets($f, 4096));
	$response['pixel'] = trim(fgets($f, 4096));
	$response['nisr'] = trim(fgets($f, 4096));
	$response['flag'] = trim(fgets($f, 4096));
	$response['machine'] = trim(fgets($f, 4096));
	$response['photons'] = trim(fgets($f, 4096));
}
else
{
	$response['ret'] = "False";
}

fclose($f);
echo json_encode($response);
?>
