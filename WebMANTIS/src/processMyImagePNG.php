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
//	Filename:	processMyImagePNG.php
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	Parses the detection data generated by VisualMANTIS and generates the PHS images
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */

$s = $_POST['n'];
$folder = $_POST['folder'];

$name = "PRF Image ".$s;
$imagename = "PRF_image".$s.".png";
$filename = $folder."/visualmantis/myimage".$s.".dat";
$f = fopen($folder."/visualmantis/myimage".$s.".dat", "r"); 

$response = array();

if($f)
{
	// writes gnuplot properties to the console
	$console = popen("gnuplot", "w");
	fwrite($console, "set title '".$name."'\n"); 
	fwrite($console, "set terminal png enh size 500, 500 crop\n");
	fwrite($console, "set out '".$imagename."'\n");
	fwrite($console, "set pm3d map\n");
	fwrite($console, "set log cb\n");
	fwrite($console, "set log z\n");
	fwrite($console, "set si sq\n");
	fwrite($console, "set cbtics format \"1e%T\"\n");
	fwrite($console, "set xlabel 'x label'\n");
	fwrite($console, "set ylabel 'y label'\n");
	fwrite($console, "unset key\n");
	fwrite($console, "set palette rgbformulae 33,13,10\n");
	fwrite($console, "splot '".$filename."'\n");
	fwrite($console, "exit\n");
	
	fclose($console);
	exec("mv ".$imagename." ".$folder."/images/");
	
    $response['ret'] = "True";
}
else
{
    $response['ret'] = "False";
}

echo json_encode($response);

?>