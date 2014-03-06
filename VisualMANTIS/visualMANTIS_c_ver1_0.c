///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
// 			     //   	       fastDETECT2 - C code  		     //
//			     //		   (optical photons transport)		     //
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
//	Associated publication: Sharma Diksha, Badal Andreu and Badano Aldo, "hybridMANTIS: a CPU-GPU Monte Carlo method for modeling indirect x-ray detectors with
//				columnar scintillators". Physics in Medicine and Biology, 57(8), pp. 2357–2372 (2012)
//
//
//	File:   	hybridMANTIS_c_ver1_0.c 			
//	Author: 	Diksha Sharma (US Food and Drug Administration)
//	Email: 		diksha.sharma@fda.hhs.gov			
//	Last updated:  	Apr 13, 2012
//
//	Modified Name:	visualMANTIS_c_ver1_0.c
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This code contains modified C code to retrieve photon information during execution of hybridMANTIS.
//			The algorithm and code is unchanged and thus is the same as the base hybridMANTIS code, however what was added was
//			code to save the Point Response Function and Pulse Height Spectrum data for each iteration. This data was parsed by
//			gnuplot scripts and displayd in the visualization.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
//
//      Header libraries
//
/////////////////////////////////////////

	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
	#include <unistd.h>
	#include <stdarg.h>

/////////////////////////////////////////
//
//      Global variables
//
/////////////////////////////////////////

	#define max_photon_per_EDE 900000	// maximum number of optical photons that can be generated per energy deposition event (EDE)

	#ifndef USING_CUDA
		#define mybufsize 2304000	// CPU buffer size: # of events sent to the CPU
	#endif

/////////////////////////////////////////
//
//      Include kernel program
//
/////////////////////////////////////////
	#include "visual_kernel_cuda_c_ver1_0.cu"

int ccounter = 1;

////////////////////////////////////////////////////////////////////////////
//				MAIN PROGRAM			          //
////////////////////////////////////////////////////////////////////////////

#ifndef USING_CUDA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// cpuoptical(): Performs optical transport using CPU 
//	  	 Input arguments: gpusize
//
// 		 gpusize - buffer size already processed in the GPU. CPU runs optical transport for rest of the buffer.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void cpuoptical_(int *gpusize)				
	{    
		printf("cpuoptical: %d\n", ccounter);
		//visual
		FILE *fp;
		char detectName[100];
		char filename[100];
		char pHistName[100];

		// command line arguments
		float xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, lbound_x, lbound_y, ubound_x, ubound_y, d_max, yield, sensorRefl;
		int pixelsize, num_primary, min_optphotons, max_optphotons, num_bins;

		float dcos[3]={0}; 		// directional cosines
		float normal[3]={0}; 		// normal to surface in case of TIR
		float pos[3] = {0}; 		// intersecting point coordinates
		float old_pos[3] = {0};  	// source coordinates
	
		int nbytes = (optical_.myctropt - (*gpusize))*sizeof(struct start_info);
		struct start_info *structa;
		structa = (struct start_info*) malloc(nbytes);
		if( structa == NULL )
			printf("\n Struct start_info array CANNOT BE ALLOCATED - %d !!", (optical_.myctropt-(*gpusize)));

		// get cpu time
		clock_t start, end;
		float num_sec;

		// get current time stamp to initialize seed input for RNG
		time_t seconds;
		seconds = time (NULL);
		struct timeval tv;
		gettimeofday(&tv,NULL);

		float rr=0.0f, theta=0.0f;	
		float r=0.0f;			// random number
		float norm=0.0f;
		int jj=0;
		int my_index=0;
		int penindex = 0;		// equal to *gpusize
		int result_algo = 0;
		unsigned long long int *num_rebound;

		// initialize random number generator (RANECU)
		int seed_input = 271828182 ; 		// ranecu seed input
		int seed[2];

		// GNU scientific library (gsl) variables
		const gsl_rng_type * Tgsl;
		gsl_rng * rgsl;
		double mu_gsl;	
				
		// output image variables
		int xdim = 0;
		int ydim = 0;
		int indexi=0, indexj=0;

		// copy to local variables from PENELOPE buffers
		xdetector = inputargs_.detx;		// x dimension of detector (in um). x in (0,xdetector)
		ydetector = inputargs_.dety;		// y dimension of detector (in um). y in (0,ydetector)
		height = inputargs_.detheight;		// height of column and thickness of detector (in um). z in range (-H/2, H/2)
		radius = inputargs_.detradius;		// radius of column (in um).
		n_C = inputargs_.detnC;			// refractive index of columns
		n_IC = inputargs_.detnIC;		// refractive index of intercolumnar material
		top_absfrac = inputargs_.dettop;	// column's top surface absorption fraction (0.0, 0.5, 0.98)
		bulk_abscoeff = inputargs_.detbulk;	// column's bulk absorption coefficient (in um^-1) (0.001, 0.1 cm^-1) 
		beta = inputargs_.detbeta;		// roughness coefficient of column walls
		d_min = inputargs_.detdmin;		// minimum distance a photon can travel when transmitted from a column
		d_max = inputargs_.detdmax;
		lbound_x = inputargs_.detlboundx;	// x lower bound of region of interest of output image (in um)
		lbound_y = inputargs_.detlboundy;	// y lower bound (in um)
		ubound_x = inputargs_.detuboundx;	// x upper bound (in um) 
		ubound_y = inputargs_.detuboundy;	// y upper bound (in um)
		yield = inputargs_.detyield;		// yield (/eV)
		pixelsize = inputargs_.detpixel;	// 1 pixel = pixelsize microns (in um)
		sensorRefl = inputargs_.detsensorRefl;	// Non-Ideal sensor reflectivity (%)
		num_primary = inputargs_.mynumhist;	// total number of primaries to be simulated
		min_optphotons = inputargs_.minphotons;	// minimum number of optical photons detected to be included in PHS
		max_optphotons = inputargs_.maxphotons;	// maximum number of optical photons detected to be included in PHS
		num_bins = inputargs_.mynumbins;	// number of bins for genrating PHS

	      	// create a generator chosen by the environment variable GSL_RNG_TYPE
	       	gsl_rng_env_setup();	     
	       	Tgsl = gsl_rng_default;
	       	rgsl = gsl_rng_alloc (Tgsl);
	       	
	       	// dimensions of PRF image
		xdim = ceil((ubound_x - lbound_x)/pixelsize);
		ydim = ceil((ubound_y - lbound_y)/pixelsize);
		unsigned long long int myimage[xdim][ydim];
		
		//initialize the output image 2D array
		for(indexi = 0; indexi < xdim; indexi++)
		{ 
		 for(indexj = 0; indexj < ydim; indexj++)
		 {
		  myimage[indexi][indexj] = 0;
		 }
		}
	
		// memory for storing histogram of # photons detected/primary
		int *h_num_detected_prim = 0;		
		h_num_detected_prim = (int*)malloc(sizeof(int)*num_primary);
			
		for(indexj=0; indexj < num_primary; indexj++)
		  h_num_detected_prim[indexj] = 0;
		  	
		// memory for storing histogram of # photons detected/primary
		int *h_histogram = 0;		
		h_histogram = (int*)malloc(sizeof(int)*num_bins);
			
		for(indexj=0; indexj < num_bins; indexj++)
		  h_histogram[indexj] = 0;

	penindex = *gpusize;
	
	for(my_index = 0; my_index < (optical_.myctropt-(*gpusize)); my_index++)		// iterate over x-rays
	{

		// reset the global counters
		num_generated = 0;
		num_detect=0;
		num_abs_top=0;	
		num_abs_bulk=0;	
		num_lost=0;
		num_outofcol=0;
		num_theta1=0;
		photon_distance=0.0f;

		//re-initialize the output image 2D array
		for(indexi = 0; indexi < xdim; indexi++)
		 for(indexj = 0; indexj < ydim; indexj++)
		   myimage[indexi][indexj] = 0;

		// copying PENELOPE buffer into *structa

		// units in the penelope output file are in cm. Convert them to microns.
		structa[my_index].str_x = optical_.xbufopt[penindex+my_index] * 10000.0f;	// x-coordinate of interaction event
		structa[my_index].str_y = optical_.ybufopt[penindex+my_index] * 10000.0f;	// y-coordinate
		structa[my_index].str_z = optical_.zbufopt[penindex+my_index] * 10000.0f;	// z-coordinate
		structa[my_index].str_E = optical_.debufopt[penindex+my_index];			// energy deposited
		structa[my_index].str_histnum = optical_.nbufopt[penindex+my_index];		// x-ray history

		// sample # optical photons based on light yield and energy deposited for this interaction event (using Poisson distribution)
		mu_gsl = (double)structa[my_index].str_E * yield;
		structa[my_index].str_N = gsl_ran_poisson(rgsl,mu_gsl);

		if(structa[my_index].str_N > max_photon_per_EDE)
		{
			printf("\n\n str_n exceeds max photons. program is exiting - %d !! \n\n",structa[my_index].str_N);
			exit(0);
		}


		num_rebound = (unsigned long long int*) malloc(structa[my_index].str_N*sizeof(unsigned long long int));
		if(num_rebound == NULL)
		printf("\n Error allocating num_rebound memory !\n");

		// start the clock
		start = clock();

		// initialize the RANECU generator in a position far away from the previous history:
		seed_input = (int)(seconds/3600+tv.tv_usec);			// seed input=seconds passed since 1970+current time in micro secs
		init_PRNG(my_index, 50000, seed_input, seed);      		// intialize RNG

		for(jj=0; jj<structa[my_index].str_N; jj++)
			num_rebound[jj] = 0;

		// reset the directional cosine and normal vectors
		dcos[0]=0.0f; dcos[1]=0.0f; dcos[2]=0.0f;
		normal[0]=0.0f; normal[1]=0.0f; normal[2]=0.0f;

		// re-initialize myimage
		for(indexi = 0; indexi < xdim; indexi++)
		 for(indexj = 0; indexj < ydim; indexj++)
			  myimage[indexi][indexj] = 0;
		
		// set starting location of photon
		pos[0] = structa[my_index].str_x; pos[1] = structa[my_index].str_y; pos[2] = structa[my_index].str_z;	
		old_pos[0] = structa[my_index].str_x; old_pos[1] = structa[my_index].str_y; old_pos[2] = structa[my_index].str_z;

		// initializing the direction cosines for the first particle in each core
		r = (ranecu(seed) * 2.0f) - 1.0f; 	// random number between (-1,1)
		 	
		while(fabs(r) <= 0.01f)	
		 {
		   	r = (ranecu(seed) * 2.0f) - 1.0f;  	
		 }

		dcos[2] = r;		// random number between (-1,1)
		rr = sqrt(1.0f-r*r);
		theta=ranecu(seed)*twopipen;
		dcos[0]=rr*cos(theta);
		dcos[1]=rr*sin(theta);

		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))	// normalize
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 }


		local_counter=0;	// total number of photons terminated (either detected at sensor, absorbed at the top or in the bulk) [global variable]
		while(local_counter < structa[my_index].str_N)		// until all the optical photons are not transported
		 { 
			
			absorbed = 0;
			detect = 0;
			bulk_abs = 0;

			// set starting location of photon
			pos[0] = structa[my_index].str_x; pos[1] = structa[my_index].str_y; pos[2] = structa[my_index].str_z;	
			old_pos[0] = structa[my_index].str_x; old_pos[1] = structa[my_index].str_y; old_pos[2] = structa[my_index].str_z;
			num_generated++;
			result_algo = 0;

			while(result_algo == 0)
			 {
			  	result_algo = algo(normal, old_pos, pos, dcos, num_rebound, seed, structa[my_index], &myimage[0][0], xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, sensorRefl, d_max, ydim, h_num_detected_prim);      
			 }

		 }	


		// end the clock		
		end = clock();

		num_sec = (((float)end - start)/CLOCKS_PER_SEC);

			
		 // type cast unsigned long long int to double
		 double cast_num_generated;
		 double cast_num_detect;
		 double cast_num_abs_top;
		 double cast_num_abs_bulk;
		 double cast_num_lost;
		 double cast_num_outofcol;
		 double cast_num_theta1;
		 double cast_gputime;

		 cast_num_generated = (double)num_generated;
		 cast_num_detect    = (double)num_detect;
		 cast_num_abs_top   = (double)num_abs_top;
		 cast_num_abs_bulk  = (double)num_abs_bulk;
		 cast_num_lost      = (double)num_lost;
		 cast_num_outofcol  = (double)num_outofcol;
		 cast_num_theta1    = (double)num_theta1;
		 cast_gputime	    = (double)(num_sec*1000.0);		// convert in millisecond

		 // save to global counters
		 optstats_.glgen      = optstats_.glgen      + cast_num_generated;
		 optstats_.gldetect   = optstats_.gldetect   + cast_num_detect;
		 optstats_.glabstop   = optstats_.glabstop   + cast_num_abs_top;
		 optstats_.glabsbulk  = optstats_.glabsbulk  + cast_num_abs_bulk;
		 optstats_.gllost     = optstats_.gllost     + cast_num_lost;
		 optstats_.gloutofcol = optstats_.gloutofcol + cast_num_outofcol;
		 optstats_.gltheta1   = optstats_.gltheta1   + cast_num_theta1;
		 optstats_.glgputime  = optstats_.glgputime  + cast_gputime;

		// release resources
		free(num_rebound);
	
	}	// my_index loop ends
	
	// make histogram of number of detected photons/primary for num_bins
	int binsize=0, newbin=0;
	int bincorr=0;
							
	binsize = floor((max_optphotons-min_optphotons)/num_bins);	// calculate size of each bin. Assuming equally spaced bins.
	bincorr = floor(min_optphotons/binsize);			// correction in bin number if min_optphotons > 0.
			
	for(indexi = 0; indexi < num_primary; indexi++)
	 {
	 	newbin = floor(h_num_detected_prim[indexi]/binsize) - bincorr;	// find bin #
		 	
	 	if(h_num_detected_prim[indexi] > 0)	// store only non-zero bins
	 	{
		 	if(h_num_detected_prim[indexi] <= min_optphotons)	// # detected < minimum photons given by user, add to the first bin
		 		h_histogram[0]++;
		 	else if(h_num_detected_prim[indexi] >= max_optphotons)	// # detected > maximum photons given by user, then add to the last bin
		 		h_histogram[num_bins-1]++;
		 	else
		 		h_histogram[newbin]++; 
		}
	 }
			
	sprintf(detectName, "detect%d.dat", ccounter);
	fp=fopen(detectName, "w");
	// add num_detected_primary to gldetprimary array in PENELOPE
	for(indexi = 0; indexi < num_bins; indexi++)
	{
		outputdetprim_.gldetprimary[indexi] = outputdetprim_.gldetprimary[indexi] + h_histogram[indexi];
		fprintf(fp, "        %d\n", outputdetprim_.gldetprimary[indexi]);
	}
	fclose(fp);

	sprintf(filename, "myimage%d.dat", ccounter);
	fp=fopen(filename, "w");
	for(indexi = 0; indexi < xdim; indexi++)
	{
		for(indexj = 0; indexj < ydim; indexj++)
		{
			fprintf(fp, "       %.4f\n",  outputimage_.newimageopt[indexi][indexj]);
		}
	}
	fclose(fp);

	// release resources
	free(structa);
	free(h_num_detected_prim);
	free(h_histogram);
	
	ccounter += 1;
	
		return;
	}	// C main() ends
	
#endif


