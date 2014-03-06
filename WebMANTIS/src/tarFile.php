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
//	Filename:	tarFile.php
//	Updated: 	2/27/2014
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This code gets all the image and data files from the various folders and zips it all together for download.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
 
$folder = $_POST['folder'];
$fname = $_POST['fname'];

//echo $file;
exec("rm -rf ".$folder."/downloads/*");
exec("cp -rf ".$folder."/images/* ".$folder."/downloads/");
exec("cp -rf ".$folder."/data/* ".$folder."/downloads/");
exec("rm -rf ".$folder."/visualmantis/detected_*");
exec("rm -rf ".$folder."/visualmantis/myimage_*");
$str = "cp -rf ".$folder."/visualmantis/*.dat ".$folder."/downloads/";
exec($str);

$phar = new PharData($fname);
// add all files in the project
$phar->buildFromDirectory($folder."/downloads/");

exec("mv ".$fname." ".$folder."/");

$response = array();
$response['filename'] = $fname;//$fname;

echo json_encode($response);

?>
