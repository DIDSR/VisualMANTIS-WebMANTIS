///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
//			     //	    	    fastDETECT2  - CUDA code 		     //
//			     //		   (optical photons transport)		     //
//			     //							     //
//			     //////////////////////////////////////////////////////////
//
// 
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
//	File:   	hybridMANTIS_cuda_ver1_0.cu 			
//	Author: 	Diksha Sharma (US Food and Drug Administration)
//	Email: 		diksha.sharma@fda.hhs.gov			
//	Last updated:  	Apr 18, 2012
// 
//	Modified Name:	visualMANTIS_cuda_ver1_0.cu
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This code contains modified C code that called CUDA kernels in order to retrieve photon information during execution of hybridMANTIS.
//			The algorithm and code is unchanged and thus is the same as the base hybridMANTIS code, however what was added was
//			additional data structures to save photon data during the execution so that it can be used by the visualization for
//			rendering. At the end of each iteration, the data structures are read back from the GPU and parsed and saved to text
//			files that are then used for visualization.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////
//
//      Header libraries
//
/////////////////////////////////////////
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef USING_CUDA
	#include <cutil_inline.h>
	#include <vector_types.h>
	#include <stdint.h>
#endif

#include <unistd.h>
#include <stdarg.h>

/////////////////////////////////////////
//
//      Global variables
//
/////////////////////////////////////////
#define max_photon_per_EDE 900000	// maximum number of optical photons that can be generated per energy deposition event (EDE)

#ifdef USING_CUDA
	#define gpubufsize 2304000	// GPU buffer size: # of events sent to the GPU
#endif

/////////////////////////////////////////
//
//      Include kernel program
//
/////////////////////////////////////////
#include "visual_kernel_cuda_c_ver1_0.cu"

/////////////////////////////////////////
//
//      CUDA parameters
//
/////////////////////////////////////////
#ifdef USING_CUDA
	#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

	#define GRIDSIZE 18000		// number of blocks
	#define BLOCKSIZE 128		// number of threads
#endif

int counter = 0;
unsigned long long int *han_h_myimage;
unsigned long long int *han_h_photonHist;
unsigned long long int *han_d_photonHist;
int * vis_num_detected_primary;
int numPhotonHist;
int finishPhotonHist = 0;
int boolToCollectData = 0;

////////////////////////////////////////////////////////////////////////////
//				MAIN PROGRAM			          //
////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// gpuoptical(): Performs optical transport using GPU 
//	  	 Input arguments: penctr, myfactGPU
//
// 		 penctr - flag to indicate how optical transport will be run on GPU. 
//		   a value of '99' : calling gpuoptical() first time to initialize the GPU and allocate memories and reset counters.
//		   a value of '100': calling optical transport kernel
//		   a value of '101': calling gpuoptical() last time; running optical transport for remaining buffer; copying data from device to host and getting output images.
// 		 myfactGPU - buffer size to be sent to GPU after load balancing
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	extern "C" void gpuoptical_(int *penctr, int *myfactGPU)	
	{    
		// command line arguments
		float xdetector, ydetector, radius, height, n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, lbound_x, lbound_y, ubound_x, ubound_y, d_max, yield, sensorRefl;
		int pixelsize, num_primary, min_optphotons, max_optphotons, num_bins;
		int i, j;
		
		// CUDA variables	
		unsigned int hTimer = 0;		// timer
		static float totalgpuTime = 0.0f;	// total time taken
		dim3 threads, blocks;			// threads and blocks
		int devID;				// GPU device ID
		
		// GNU scientific library (gsl) variables
		const gsl_rng_type * Tgsl;
		gsl_rng * rgsl;
		double mu_gsl;	
		
		// Host (CPU) counters		
		unsigned long long int host_num_generated = 0; 	// total # of photons generated for all the x-ray histories (across all threads)
		unsigned long long int host_num_detect = 0;	// total # of photons detected at the sensor plane of a column
		unsigned long long int host_num_lost = 0;	// total # of photons lost when exiting out of the detector boundaries in x/y direction
		unsigned long long int host_num_abs_top = 0;	// total # of photons absorbed at the top surface of the detector
		unsigned long long int host_num_abs_bulk = 0;	// total # of photons absorbed in the bulk of the detector
		unsigned long long int host_num_outofcol = 0;	// total # of photons killed because they moved out of current column when reflected (due to precision errors)
		unsigned long long int host_num_theta1 = 0;	// total # of photons killed if incidence angle > 1.57 or < 0 radian (after resampling max 100 times)	

		// create 2D array for storing output PRF image
		int xdim = 0;
		int ydim = 0;
		int indexi=0, indexj=0;
		int my_index=0;
		size_t pitch;				// pitch used for storing 2D image array in GPU memory

		int nbytes = (*myfactGPU)*sizeof(struct start_info);	// total number of bytes for storing interaction events buffer information

		// allocate memory pointers
		unsigned long long int *myimage = 0;	// device memory for output image
		int *num_detected_primary = 0;		// device memory for # detected photons/primary
		struct start_info *h_a = 0;             // pointer to the struct info data in the host memory
		struct start_info *d_a = 0;             // pointers to struct data in the device memory

		// copy to local variables from PENELOPE common block
		xdetector = inputargs_.detx;		// x dimension of detector (in um). x in (0,xdetector)
		ydetector = inputargs_.dety;		// y dimension of detector (in um). y in (0,ydetector)
		height = inputargs_.detheight;		// height of column and thickness of detector (in um). z in range (-H/2, H/2)
		radius = inputargs_.detradius;		// radius of column (in um).
		n_C = inputargs_.detnC;			// refractive index of columns
		n_IC = inputargs_.detnIC;		// refractive index of intercolumnar material
		top_absfrac = inputargs_.dettop;	// column's top surface absorption fraction (0,1)
		bulk_abscoeff = inputargs_.detbulk;	// column's bulk absorption coefficient (in um^-1) 
		beta = inputargs_.detbeta;		// roughness coefficient of column walls (0,0.5)
		d_min = inputargs_.detdmin;		// minimum distance a photon can travel when transmitted from a column
		d_max = inputargs_.detdmax;		// maximum distance a photon can travel when transmitted from a column
		lbound_x = inputargs_.detlboundx;	// x lower bound of region of interest of output PRF image (in um)
		lbound_y = inputargs_.detlboundy;	// y lower bound (in um)
		ubound_x = inputargs_.detuboundx;	// x upper bound (in um) 
		ubound_y = inputargs_.detuboundy;	// y upper bound (in um)
		yield = inputargs_.detyield;		// light yield (/eV)
		pixelsize = inputargs_.detpixel;	// 1 pixel = 'pixelsize' microns (in um)
		sensorRefl = inputargs_.detsensorRefl;	// Non-Ideal sensor reflectivity (%) (0,100)
		num_primary = inputargs_.mynumhist;	// total number of primaries to be simulated
		min_optphotons = inputargs_.minphotons;	// minimum number of optical photons detected to be included in PHS
		max_optphotons = inputargs_.maxphotons;	// maximum number of optical photons detected to be included in PHS
		num_bins = inputargs_.mynumbins;	// number of bins for genrating PHS
				
		// set the device with max GFlops	
		devID = cutGetMaxGflopsDeviceId();
		cudaSetDevice( devID );

	    // create a generator chosen by the 
		//  environment variable GSL_RNG_TYPE 
	    gsl_rng_env_setup();	     
	    Tgsl = gsl_rng_default;
	    rgsl = gsl_rng_alloc (Tgsl);

		// dimensions of PRF image
		xdim = ceil((ubound_x - lbound_x)/pixelsize);
		ydim = ceil((ubound_y - lbound_y)/pixelsize);

		//han
		unsigned long long int *h_myimage = 0;          // page-locked host memory for asynchronous copying (contain host image for evey kernel run)
		int *h2_num_detected_primary = 0;
		struct histStruct * hPhotonHist;
		struct histStruct * dPhotonHist;
		printf("%d\n", *penctr);

		if(*penctr == 99)	// initialize GPU; allocate and initialize memories
		{
			// allocate device memory for storing output arrays
			cudaMallocPitch((void**)&myimage, &pitch, xdim*sizeof(unsigned long long int), ydim);		// allocate 2D image array (PRF)
			cutilSafeCall( cudaMemset2D(myimage, pitch, 0, xdim*sizeof(unsigned long long int), ydim) );	// initialize to 0

			cutilSafeCall( cudaMalloc((void**)&num_detected_primary, sizeof(int)*num_primary) );		// create 1D array for outputting # detected/primary
			cutilSafeCall( cudaMemset(num_detected_primary, 0, sizeof(int)*num_primary) );			// initialize to 0

			// allocate host and device memory for transferring buffer information
			cutilSafeCall( cudaMallocHost((void**)&h_a, nbytes) ); 
			cutilSafeCall( cudaMalloc((void**)&d_a, nbytes) );

			// copy address of memory pointers to PENELOPE variables. These variables are used later to point to device and host memories without re-initializing the GPU.
			gpumemaddr_.gpuimage 	= (unsigned long long int)myimage;
			gpumemaddr_.gpudetect 	= (unsigned long long int)num_detected_primary;
			gpumemaddr_.hosta 	= (unsigned long long int)h_a;	
			gpumemaddr_.deva 	= (unsigned long long int)d_a;
		    gpumemaddr_.devpitch 	= (unsigned long long int)pitch;

			// reset the host counters
			host_num_generated=0;
			host_num_detect=0;
			host_num_abs_top=0;	
			host_num_abs_bulk=0;	
			host_num_lost=0;
			host_num_outofcol=0;
			host_num_theta1=0;

			FILE * fphist;
			char numPhotonHistStr[100];
			fphist = fopen("numPhotonHistories.txt", "r");
			fgets(numPhotonHistStr, 100, fphist);
			numPhotonHist = atoi(numPhotonHistStr);
			fclose(fphist);
			
			// reset device counters to zero
			cutilSafeCall(cudaMemcpyToSymbol("num_detect",&host_num_detect,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
			cutilSafeCall(cudaMemcpyToSymbol("num_generated",&host_num_generated,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
			cutilSafeCall(cudaMemcpyToSymbol("num_abs_top",&host_num_abs_top,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
			cutilSafeCall(cudaMemcpyToSymbol("num_abs_bulk",&host_num_abs_bulk,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));	
			cutilSafeCall(cudaMemcpyToSymbol("num_lost",&host_num_lost,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));
			cutilSafeCall(cudaMemcpyToSymbol("num_outofcol",&host_num_outofcol,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));
			cutilSafeCall(cudaMemcpyToSymbol("num_theta1",&host_num_theta1,sizeof(unsigned long long int)*1,0,cudaMemcpyHostToDevice));
			cutilSafeCall(cudaMemcpyToSymbol("dev_numPhotonHist", &numPhotonHist, sizeof(int)*1, 0 , cudaMemcpyHostToDevice));

			//han
			//allocate memory for h_myimage
			cutilSafeCall( cudaMallocHost((void**)&h_myimage, xdim*ydim*sizeof(unsigned long long int)) ); 
			//allocate memory for h2_num_detected_primary
			cutilSafeCall( cudaMallocHost((void**)&h2_num_detected_primary, sizeof(int) * num_primary) ); 
			
			// allocate host and device memory for transferring buffer information in photon histories
			cutilSafeCall( cudaMallocHost((void**)&hPhotonHist, numPhotonHist*sizeof(struct histStruct)) ); 
			cutilSafeCall( cudaMalloc((void**)&dPhotonHist, numPhotonHist*sizeof(struct histStruct)) );
			
			for(i=0;i<numPhotonHist;i++)
			{
				hPhotonHist[i].histCounter = 0;
			}
			cutilSafeCall( cudaMemcpy(dPhotonHist, hPhotonHist, numPhotonHist*sizeof(struct histStruct), cudaMemcpyHostToDevice) );
			
			// copy memory pointers to avoid re mallocing
			han_h_myimage = h_myimage;
			han_h_photonHist = (unsigned long long int *)hPhotonHist;
			han_d_photonHist = (unsigned long long int *)dPhotonHist;
			vis_num_detected_primary = h2_num_detected_primary;
			
			cutilCheckError( cutCreateTimer(&hTimer) );

		}
		else if(*penctr == 100)			// run optical kernel
		{
			FILE *fp;
			char detectName[100];
			char filename[100];
			char pHistName[100];
				
			// synchronize threads to ensure that GPU is not busy processing previous kernel call
			cudaThreadSynchronize();
			hPhotonHist = (struct histStruct *) han_h_photonHist;
			dPhotonHist = (struct histStruct *) han_d_photonHist;
			
			if(counter == 0)
			{
				counter +=1;
			}
			else
			{
				h_myimage = han_h_myimage;
				h2_num_detected_primary = vis_num_detected_primary;
				
				sprintf(filename, "myimage%d.dat", counter);
				fp=fopen(filename, "w");
				
				// add h_myimage to the new_myimage (array in PENELOPE)
				for(indexi = 0; indexi < ydim; indexi++)
				{
		 	 		for(indexj = 0; indexj < xdim; indexj++)
					{
						outputimage_.newimageopt[indexi][indexj] = outputimage_.newimageopt[indexi][indexj] + h_myimage[indexi*xdim + indexj];
						fprintf(fp, "       %.4f\n",  outputimage_.newimageopt[indexi][indexj] * (1.0 / optical_.nbufopt[(*myfactGPU)-1]));
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
				
				/*for(i=0;i<numPhotonHist;i++)
				{
					sprintf(pHistName, "photon_hist_%d_%d.dat", counter, i);
				
					fp = fopen(pHistName, "w");	
					fprintf(fp, "%d\n", hPhotonHist[i].histCounter);
					for(j=0;j<hPhotonHist[i].histCounter;j++)
					{
						fprintf(fp, "%f %f %f %f %f %d\n", hPhotonHist[i].x[j], hPhotonHist[i].y[j], hPhotonHist[i].z[j], hPhotonHist[i].Xc[j], hPhotonHist[i].Yc[j], hPhotonHist[i].terminated[j]);
					}
					
					fclose(fp);
				}*/
				
				int indexi, indexj;
				int *h_histogram = 0;		// host memory for storing histogram of # photons detected/primary
				h_histogram = (int*)malloc(sizeof(int)*num_bins);
				
				for(indexj=0; indexj < num_bins; indexj++)
				{
					h_histogram[indexj] = 0;
				}
				
				// make histogram of number of detected photons/primary for num_bins
				int binsize=0, newbin=0;
				int bincorr=0;
							
				binsize = floor((max_optphotons-min_optphotons)/num_bins);	// calculate size of each bin. Assuming equally spaced bins.
				bincorr = floor(min_optphotons/binsize);			// correction in bin number if min_optphotons > 0.
			
				for(indexi = 0; indexi < num_primary; indexi++)
				{
					newbin = floor(h2_num_detected_primary[indexi]/binsize) - bincorr;	// find bin #
			 	
					if(h2_num_detected_primary[indexi] > 0)	// store only non-zero bins
					{
						if(h2_num_detected_primary[indexi] <= min_optphotons)	// # detected < minimum photons given by user, add to the first bin
							h_histogram[0]++;
						else if(h2_num_detected_primary[indexi] >= max_optphotons)	// # detected > maximum photons given by user, then add to the last bin
							h_histogram[num_bins-1]++;
						else
							h_histogram[newbin]++; 
					}
				}
				
				sprintf(detectName, "detect%d.dat", counter);
				fp=fopen(detectName, "w");
				
				// add num_detected_primary to gldetprimary array in PENELOPE
				for(indexi = 0; indexi < num_bins; indexi++)
				{
					outputdetprim_.gldetprimary[indexi] = outputdetprim_.gldetprimary[indexi] + h_histogram[indexi];
					fprintf(fp, "        %d\n", outputdetprim_.gldetprimary[indexi]);
				}
				fclose(fp);
				
				counter += 1;
				printf("Finished writing photon history and image data.\n");
			}
			
			if(finishPhotonHist < numPhotonHist)
			{
				for(i=0;i<numPhotonHist;i++)
				{
					if(hPhotonHist[i].histCounter > 0 && finishPhotonHist < numPhotonHist)
					{
						sprintf(pHistName, "photon_hist_1_%d.dat", finishPhotonHist);
				
						fp = fopen(pHistName, "w");	
						fprintf(fp, "%d\n", hPhotonHist[i].histCounter);
						for(j=0;j<hPhotonHist[i].histCounter;j++)
						{
							fprintf(fp, "%f %f %f %f %f %d\n", hPhotonHist[i].x[j], hPhotonHist[i].y[j], hPhotonHist[i].z[j], hPhotonHist[i].Xc[j], hPhotonHist[i].Yc[j], hPhotonHist[i].terminated[j]);
						}
					
						finishPhotonHist += 1;
						fclose(fp);
					}
				}
			}
			
			// allocate nbytes
			nbytes = (*myfactGPU)*sizeof(struct start_info);

			// copy memory address from PENELOPE variables
			myimage = (unsigned long long int*)gpumemaddr_.gpuimage;
			num_detected_primary = (int*)gpumemaddr_.gpudetect;
			h_a = (struct start_info*)gpumemaddr_.hosta;
			d_a = (struct start_info*)gpumemaddr_.deva;
			pitch = (size_t)gpumemaddr_.devpitch;
			
			// assign number of threads and blocks
			threads = dim3(BLOCKSIZE,1);
			blocks = dim3(GRIDSIZE,1);

			// reading data from buffer
			for(my_index = 0; my_index < (*myfactGPU); my_index++)		// iterate over buffer length
			{
				// units in the penelope output file are in cm. Convert them to microns.
				h_a[my_index].str_x = optical_.xbufopt[my_index] * 10000.0f;	// x-coordinate of interaction event
				h_a[my_index].str_y = optical_.ybufopt[my_index] * 10000.0f;	// y-coordinate
				h_a[my_index].str_z = optical_.zbufopt[my_index] * 10000.0f;	// z-coordinate
				h_a[my_index].str_E = optical_.debufopt[my_index];		// energy deposited
				h_a[my_index].str_histnum = optical_.nbufopt[my_index];		// x-ray history number
				
				// sample # optical photons based on light yield and energy deposited for this interaction event (using Poisson distribution)
				mu_gsl = (double)h_a[my_index].str_E * yield;
				h_a[my_index].str_N = gsl_ran_poisson(rgsl,mu_gsl);
				
				if(h_a[my_index].str_N > max_photon_per_EDE)
				{
					printf("\n\n str_n exceeds max photons. program is exiting - %d !! \n\n", h_a[my_index].str_N);
					exit(0);
				}
			} 
			
			cudaMemset(dPhotonHist, 0, numPhotonHist*sizeof(struct histStruct));
			memset(hPhotonHist, 0, numPhotonHist*sizeof(struct histStruct));
			
			// execute the optical transport kernel 
			// asynchronously copy data from host to device	(all to stream 0)
			cutilSafeCall( cudaMemcpyAsync(d_a, h_a, nbytes, cudaMemcpyHostToDevice, 0) );
			
			if(finishPhotonHist < numPhotonHist)
			{
				boolToCollectData = 1;
			}
			else
			{
				boolToCollectData = 0;
			}
			
			// each kernel has BLOCKSIZE threads; each thread transports one event in the buffer (info.str_N optical photons)
			algo<<<blocks, threads, 0, 0>>>(d_a, myimage, num_detected_primary, pitch, (*myfactGPU), xdetector, ydetector, radius, height, 
			n_C, n_IC, top_absfrac, bulk_abscoeff, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, d_max, sensorRefl, dPhotonHist, boolToCollectData); 

			// asynchronously copy data from device to host
			h_myimage = han_h_myimage;
			cutilSafeCall( cudaMemcpy2DAsync((void*)h_myimage,sizeof(unsigned long long int)*xdim,(void*)myimage,pitch, sizeof(unsigned long long int)*xdim,ydim,cudaMemcpyDeviceToHost, 0) );
			han_h_myimage = h_myimage;
			
			cutilSafeCall( cudaMemcpyAsync(hPhotonHist, dPhotonHist, numPhotonHist*sizeof(struct histStruct), cudaMemcpyDeviceToHost) );
			han_h_photonHist = (unsigned long long int *)hPhotonHist;
			han_d_photonHist = (unsigned long long int *)dPhotonHist;
			
			h2_num_detected_primary = vis_num_detected_primary;
			cutilSafeCall( cudaMemcpyAsync(h2_num_detected_primary, num_detected_primary, sizeof(int)*num_primary, cudaMemcpyDeviceToHost, 0) );
			vis_num_detected_primary = h2_num_detected_primary;
			
			cutilCheckMsg("algo() execution failed\n");
		}	
		else if(*penctr == 101)		// calling optical transport kernel last time. copy back final data from device to host
		{

			// synchronize threads to ensure that GPU is not busy processing previous kernel call
			cudaThreadSynchronize();
			// here 'nbytes' is not necesaarily equal to 'gpubufsize*sizeof(struct start_info)', because in the last call optical_.myctropt can be <= gpubufsize
			nbytes = (*myfactGPU)*sizeof(struct start_info);

			// copy memory address from PENELOPE variables
			myimage = (unsigned long long int*)gpumemaddr_.gpuimage;
			num_detected_primary = (int*)gpumemaddr_.gpudetect;
			h_a = (struct start_info*)gpumemaddr_.hosta;
			d_a = (struct start_info*)gpumemaddr_.deva;
			pitch = (size_t)gpumemaddr_.devpitch;

			// allocate host memory
			//unsigned long long int *h_myimage = 0;          // page-locked host memory for asynchronous copying (contain host image for evey kernel run)
			//allocate memory for h_myimage
			//cutilSafeCall( cudaMallocHost((void**)&h_myimage, xdim*ydim*sizeof(unsigned long long int)) ); 
			h_myimage = han_h_myimage;
			hPhotonHist = (struct histStruct *) han_h_photonHist;
			dPhotonHist = (struct histStruct *) han_d_photonHist;
			
			int *h_num_detected_primary = 0;		// host memory to get # detected/primary
			cutilSafeCall( cudaMallocHost((void**)&h_num_detected_primary, sizeof(int)*num_primary) );

			for(indexj=0; indexj < num_primary; indexj++)
			  h_num_detected_primary[indexj] = 0;
			  
			int *h_histogram = 0;		// host memory for storing histogram of # photons detected/primary
			h_histogram = (int*)malloc(sizeof(int)*num_bins);
			
			for(indexj=0; indexj < num_bins; indexj++)
			  h_histogram[indexj] = 0;

			// assign number of threads and blocks
			threads = dim3(BLOCKSIZE,1);
			blocks = dim3(GRIDSIZE,1);


			// reading data from buffer
			for(my_index = 0; my_index < (*myfactGPU); my_index++)		// iterate over x-rays
			{
				// units in the penelope output file are in cm. Convert them to microns.
				h_a[my_index].str_x = optical_.xbufopt[my_index] * 10000.0f;	// x-coordinate
				h_a[my_index].str_y = optical_.ybufopt[my_index] * 10000.0f;	// y-coordinate
				h_a[my_index].str_z = optical_.zbufopt[my_index] * 10000.0f;	// z-coordinate
				h_a[my_index].str_E = optical_.debufopt[my_index];		// energy deposited
				h_a[my_index].str_histnum = optical_.nbufopt[my_index];		// x-ray history number

				// sample # optical photons based on light yield and energy deposited for this interaction event
				mu_gsl = (double)h_a[my_index].str_E * yield;
				h_a[my_index].str_N = gsl_ran_poisson(rgsl,mu_gsl);

				if(h_a[my_index].str_N > max_photon_per_EDE)
				{
					printf("\n\n str_n exceeds max photons. program is exiting - %d !! \n\n", h_a[my_index].str_N);
					exit(0);
				}

			} // for loop ends

			// execute the kernel 
			// asynchronously copy data from host to device	(all to stream 0)
			cutilSafeCall( cudaMemcpyAsync(d_a, h_a, nbytes, cudaMemcpyHostToDevice, 0) );

			// each kernel has BLOCKSIZE threads; each thread transports one event in the buffer (info.str_N optical photons)
			algo<<<blocks, threads, 0, 0>>>(d_a, myimage, num_detected_primary, pitch, (*myfactGPU), xdetector, ydetector, radius, height, n_C, n_IC, 
			top_absfrac, bulk_abscoeff, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, d_max, sensorRefl, dPhotonHist, 0); 

			// asynchronously copy image data from device to host
			cutilSafeCall( cudaMemcpy2DAsync((void*)h_myimage,sizeof(unsigned long long int)*xdim,(void*)myimage,pitch,sizeof(unsigned long long int)*xdim,ydim,cudaMemcpyDeviceToHost, 0) );
			cutilSafeCall( cudaMemcpyAsync(h_num_detected_primary, num_detected_primary, sizeof(int)*num_primary, cudaMemcpyDeviceToHost, 0) );

			cutilCheckMsg("algo() execution failed\n");

			cudaThreadSynchronize();	// ensure that GPU has finished before copying back the final results.

			// copy counters from device to host
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_detect,num_detect,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_generated,num_generated,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_abs_top,num_abs_top,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));	
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_abs_bulk,num_abs_bulk,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_lost,num_lost,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_outofcol,num_outofcol,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));
			cutilSafeCall(cudaMemcpyFromSymbol((void *) &host_num_theta1,num_theta1,sizeof(unsigned long long int)*1,0,cudaMemcpyDeviceToHost));

			// add h_myimage to the new_myimage (array in PENELOPE)
			for(indexi = 0; indexi < ydim; indexi++)
		 	 for(indexj = 0; indexj < xdim; indexj++)
				outputimage_.newimageopt[indexi][indexj] = outputimage_.newimageopt[indexi][indexj] + h_myimage[indexi*xdim + indexj];

			// make histogram of number of detected photons/primary for num_bins
			int binsize=0, newbin=0;
			int bincorr=0;
							
			binsize = floor((max_optphotons-min_optphotons)/num_bins);	// calculate size of each bin. Assuming equally spaced bins.
			bincorr = floor(min_optphotons/binsize);			// correction in bin number if min_optphotons > 0.
			
			for(indexi = 0; indexi < num_primary; indexi++)
			 {
			 	newbin = floor(h_num_detected_primary[indexi]/binsize) - bincorr;	// find bin #
			 	
 			 	if(h_num_detected_primary[indexi] > 0)	// store only non-zero bins
 			 	{
				 	if(h_num_detected_primary[indexi] <= min_optphotons)	// # detected < minimum photons given by user, add to the first bin
						h_histogram[0]++;
				 	else if(h_num_detected_primary[indexi] >= max_optphotons)	// # detected > maximum photons given by user, then add to the last bin
			 			h_histogram[num_bins-1]++;
			 		else
				 		h_histogram[newbin]++; 
				}
			 }
			
			// add num_detected_primary to gldetprimary array in PENELOPE
			for(indexi = 0; indexi < num_bins; indexi++)
				outputdetprim_.gldetprimary[indexi] = outputdetprim_.gldetprimary[indexi] + h_histogram[indexi];

		   
			// type cast unsigned long long int to double
			double cast_host_num_generated;
			double cast_host_num_detect;
			double cast_host_num_abs_top;
			double cast_host_num_abs_bulk;
			double cast_host_num_lost;
			double cast_host_num_outofcol;
			double cast_host_num_theta1;
			double cast_gputime;

			cast_host_num_generated = (double)host_num_generated;
			cast_host_num_detect    = (double)host_num_detect;
			cast_host_num_abs_top   = (double)host_num_abs_top;
			cast_host_num_abs_bulk  = (double)host_num_abs_bulk;
			cast_host_num_lost      = (double)host_num_lost;
			cast_host_num_outofcol  = (double)host_num_outofcol;
			cast_host_num_theta1    = (double)host_num_theta1;
			cast_gputime		= (double)totalgpuTime;

			 // save to global counters
			 optstats_.glgen      = optstats_.glgen      + cast_host_num_generated;
			 optstats_.gldetect   = optstats_.gldetect   + cast_host_num_detect;
			 optstats_.glabstop   = optstats_.glabstop   + cast_host_num_abs_top;
			 optstats_.glabsbulk  = optstats_.glabsbulk  + cast_host_num_abs_bulk;
			 optstats_.gllost     = optstats_.gllost     + cast_host_num_lost;
			 optstats_.gloutofcol = optstats_.gloutofcol + cast_host_num_outofcol;
			 optstats_.gltheta1   = optstats_.gltheta1   + cast_host_num_theta1;
			 optstats_.glgputime  = optstats_.glgputime  + cast_gputime;


	 
			// release resources
			cutilSafeCall(cudaFree(d_a));
			cutilSafeCall(cudaFree(myimage));
			cutilSafeCall(cudaFree(num_detected_primary));
			cutilSafeCall(cudaFree(dPhotonHist));
			cudaFreeHost(h_a);
			cudaFreeHost(h_myimage);
			cudaFreeHost(h_num_detected_primary);
			cudaFreeHost(hPhotonHist);

			free(h_histogram);
		}	// else ends


		return;
	}	// CUDA main() ends
	
#endif


