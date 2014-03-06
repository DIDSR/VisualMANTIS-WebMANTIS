///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 			     //////////////////////////////////////////////////////////
//  			     //							     //
// 			     //   	        hybridMANTIS v1.0		     //
// 			     //   	  fastDETECT2 kernel - CUDA + C  	     //
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
//	File:   	kernel_cuda_c_ver1_0.cu 			
//	Author: 	Diksha Sharma (US Food and Drug Administration)
//	Email: 		diksha.sharma@fda.hhs.gov			
//	Last updated:  	Apr 13, 2012
//
//	Modified Name:	visual_kernel_cuda_c_ver1_0.cu
//	Updated: 	4/23/2013
//	Author:		Han Dong (US Food and Drug Administration / University of Maryland Baltimore County
//	Email:		han6@umbc.edu
//	Comments:	This code contains modified CUDA kernels in order to retrieve photon information during execution of hybridMANTIS.
//			The algorithm and code is unchanged and thus is the same as the base hybridMANTIS code, however what was added was
//			additional data structures to save photon data during the execution so that it can be used by the visualization for
//			rendering.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
//
//      Header libraries
//
/////////////////////////////////////////

	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <sys/time.h>
	#include <time.h>

/////////////////////////////////////////
//
//       Constants
//
/////////////////////////////////////////

	#define twopipen 6.283185308	// 2*PI
	#define pi 3.14159265		// PI
	#define epsilon 8.1929093e-6	// a very small number for float comparisons


/////////////////////////////////////////////////////////////////////////////////////
//
//     Data structure for storing a scintillation event location and deposited energy
//
/////////////////////////////////////////////////////////////////////////////////////

	struct start_info
	{
		double str_x;		// x-coordinate
		double str_y;		// y-coordinate
		double str_z;		// z-coordinate
		double str_E;		// deposited energy
		int str_histnum;	// x-ray history #
		int str_N;		// # of optical photons to be transported for this energy
	};


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Fortran structure declarations - using PENELOPE 2006 (coded in Fortran)
//	A 'common' block in Fortran needs to be declared here to allow calling function interexchangebly.	
// 	
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Similar structure to 'start_info' - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double xbufopt[gpubufsize];	// x-coordinate array of scintillation events
			double ybufopt[gpubufsize];	// y-coordinate
			double zbufopt[gpubufsize];	// z-coordinate
			double debufopt[gpubufsize];	// deposited energy
			int nbufopt[gpubufsize];	// x-ray history #
			int myctropt;			// equal to 'gpubufsize' (buffer size when calling fastDETECT2 in the GPU)
		        int cpu_num_real;		// equal to number of x-ray histories to be run by cpuoptical
		} optical_;
	#else
		extern struct
		{
			double xbufopt[mybufsize];
			double ybufopt[mybufsize];
			double zbufopt[mybufsize];
			double debufopt[mybufsize];	
			int nbufopt[mybufsize];		
			int myctropt;			// equal to 'mybufsize'	(buffer size when calling fastDETECT2 in the CPU)
		        int cpu_num_real;
		} optical_;
	#endif

// Storing optical output statistics - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double glgen;			// total # optical photons generated
			double gldetect;		// total # optical photons detected
			double glabstop;		// total # optical photons absorbed at top surface (H/2)
			double glabsbulk;		// total # optical photons absorbed in the bulk
			double gllost;			// total # optical photons lost at detector boundaries
			double gloutofcol;		// total # optical photons killed when they went out of the current column (due to precision errors)
			double gltheta1;		// total # optical photons killed when theta1 > 90degrees; 100x resampling
			double glgputime;		// total GPU time (ms)
		} optstats_;
	#else
		extern struct
		{
			double glgen;
			double gldetect;
			double glabstop;
			double glabsbulk;
			double gllost;
			double gloutofcol;
			double gltheta1;
			double glgputime;		// total CPU time (ms)
		} optstats_;
	#endif

// Storing deposited energy and # optical photons detected - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			int gldetprimary[1000];		// total # optical photons detected per primary
		} outputdetprim_;
	#else
		extern struct
		{
			int gldetprimary[1000];
		} outputdetprim_;
	#endif

// Structure for storing point response functions - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			unsigned long long int newimageopt[501][501];	// storing point response function
			unsigned long long int tempimageopt[501][501];	// used in timing CPU optical
		} outputimage_;
	#else
		extern struct
		{
			unsigned long long int newimageopt[501][501];
			unsigned long long int tempimageopt[501][501];
		} outputimage_;
	#endif

// Storing the memory addresses of arrays, in order to call fastDETECT2 in the GPU asynchronously - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			unsigned long long int gpuimage;		// storing memory address of gpu image array (2D)
			unsigned long long int gpudetect;		// storing memory address of gpu detected array
			unsigned long long int hosta;			// storing memory address of host memory for x,y,z,E
			unsigned long long int deva;			// storing memory address of device memory for x,y,z,E
			unsigned long long int devpitch;		// storing memory address pitch of gpu image array (2D)
		} gpumemaddr_;
	#endif

// Storing the input arguments - declared in PENELOPE
	#ifdef USING_CUDA
		extern "C" struct
		{
			double detx;		// detector length in x (in um). x in (0,detx).
			double dety;		// detector length in y (in um). y in (0,dety).
			double detheight;	// height of a column or thickness of detector (in um). z in range (-H/2, H/2).
			double detradius;	// radius of a column (in um). Assuming same properties for all the columns.
			double detnC;		// refractive index of columns.
			double detnIC;		// refractive index of intercolumnar material.
			double dettop;		// column's top surface absorption fraction (0,1).
			double detbulk;		// column's bulk absorption coefficient (in um^-1). 
			double detbeta;		// roughness coefficient of column surface walls (0,0.5).
			double detdmin;		// minimum distance a photon can travel when transmitted from a column: on-the-fly geometry.
			double detdmax;		// maximum distance a photon can travel when transmitted from a column: on-the-fly geometry.
			double detlboundx;	// x lower bound of point response function (in um).
			double detlboundy;	// y lower bound (in um).
			double detuboundx;	// x upper bound (in um). 
			double detuboundy;	// y upper bound (in um). 
			double detyield;	// light yield (/eV).
			double detsensorRefl;	// Non-Ideal sensor reflectivity (%) (0,100).
			int detpixel;		// 1 pixel = 'pixelsize' microns (in um) - for storing PRF.
			int rungpu;		// flag to run on the GPU (value=1).
			int machinenum;		// machine number where code executes. This number is useful in differentiating output file names when using same input arguments.
			int mynumhist;		// total number of primaries to be simulated
			int minphotons;		// minimum number of optical photons detected to be included in PHS
			int maxphotons;		// maximum number of optical photons detected to be included in PHS
			int mynumbins;		// number of bins for genrating PHS
		} inputargs_;
	#else
		extern struct
		{
			double detx;
			double dety;
			double detheight;
			double detradius;
			double detnC;
			double detnIC;
			double dettop;
			double detbulk;
			double detbeta;
			double detdmin;
			double detdmax;
			double detlboundx;
			double detlboundy;
			double detuboundx;
			double detuboundy;
			double detyield;
			double detsensorRefl;
			int detpixel;
			int rungpu;
			int machinenum;
			int mynumhist;
			int minphotons;		
			int maxphotons;		
			int mynumbins;		
		} inputargs_;
	#endif

/////////////////////////////////////////////////////////////////////////////////////
//
//     Data structure for storing photon histories
//
/////////////////////////////////////////////////////////////////////////////////////
struct histStruct
{
	float x[1000];			// x-coordinate
	float y[1000];			// y-coordinate
	float z[1000];			// z-coordinate
	float Xc[1000];			// x coordinate of cylinder
	float Yc[1000];			// y coordinate of cylinder
	int terminated[1000];	//counter to indicate type of termination
	int histCounter;	// counter to keep track of number of histories
};
	
/////////////////////////////////////////
//
//       Function declarations
//
/////////////////////////////////////////

// transports optical photon from its generation until it ends (detected/absorbed/lost).
	#ifdef USING_CUDA
		__global__ void algo(struct start_info *info, unsigned long long int *myimage, int *num_detected_primary, size_t pitch, int rowsread, float xdetector, float ydetector, 
		float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, 
		float ubound_y, float d_max, float sensorRefl, struct histStruct * dPhotonHist, int boolToCollectData);	
	#else
		int algo(float *normal, float *old_pos, float *pos, float *dcos, unsigned long long int *num_rebound, int* seed, struct start_info info, 
		unsigned long long int *myimage, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, 
		float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, float sensorRefl, float d_max, int ydim, int *h_num_detected_prim);
	#endif

// photon within a column. calculate if it gets absorbed or moves inside the column.
	#ifdef USING_CUDA
		__device__ inline int isotropic(float3 *pos, float3 *dcos, int2* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector, 
		struct start_info *info, unsigned long long int mynum_rebound, float *Xc, float *Yc, int mytid, struct histStruct * dPhotonHist, int boolToCollectData);
	#else
		int isotropic(float *pos, float *dcos, int* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector, struct start_info info, 
		unsigned long long int mynum_rebound);
	#endif

// photon within a column. calculate distance to next position in the same column and move it.
	#ifdef USING_CUDA
		__device__ float dist_to_surface(float3 *pos, float3 *dcos, float R, float H, float xdetector, float ydetector, struct start_info *info, 
		unsigned long long int mynum_rebound, float *Xc, float *Yc, int mytid, struct histStruct * dPhotonHist, int boolToCollectData);	
	#else
		float dist_to_surface(float *pos, float *dcos, float R, float H, float xdetector, float ydetector, struct start_info info, unsigned long long int mynum_rebound);
	#endif

// photon within/between columns. calculate if it gets reflected or transmitted.
	#ifdef USING_CUDA
		__device__ int boundary_analysis(float3 *normal, float3 *pos, float3 *dcos, int2* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, 
		float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, 
		float *Xc, float *Yc, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float d_max, float sensorRefl, struct histStruct * dPhotonHist,
		int boolToCollectData);	
	#else
		int boundary_analysis(float *normal, float *pos, float *dcos, int* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, 
		float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, 
		float d_max, float sensorRefl, int ydim, int *h_num_detected_prim);	
	#endif

// transmit photon to another column. calculates the new position when it transmits, build new column and move photon here.
	#ifdef USING_CUDA
		__device__ int transmit(float3 *pos, float3 *dcos, float3 *normal, int2* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, int pixelsize, 
		float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, 
		int *num_detected_primary, float d_max, float sensorRefl, int flagCCT, float *Xc, float *Yc, struct histStruct * dPhotonHist, int boolToCollectData);
	#else
		int transmit(float *pos, float *dcos, float *normal, int* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x,
		float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, float sensorRefl, int ydim, 
		int flagCCT, int *h_num_detected_prim);	
	#endif

// called when photon reflects from sensor plane (bottom surface) of the detector, outside of any column.
	#ifdef USING_CUDA
		__device__ int refl_bottom(float3 *pos, float3 *dcos, float3 *normal, float xdetector, float ydetector, int2* seed, 
		float beta, float d_min, float H, float d_max, int mytid, float *Xc, float *Yc, struct histStruct * dPhotonHist, int boolToCollectData);
	#else
		int refl_bottom(float *pos, float *dcos, float *normal, float xdetector, float ydetector, int* seed, float beta, float d_min, float H, float d_max);	
	#endif

// calculate dot product of two vectors to give cosine of angle between them.
	#ifdef USING_CUDA
		__device__ inline float dot_product(float3 *aa, float3 *b);	
	#else
		float dot_product(float *aa, float *b);	
	#endif

// determine if photon got detected at sensor plane or is reflected back within the column
	#ifdef USING_CUDA
		__device__ inline int detection(float3 *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float sensorRefl, float d_min, int2* seed, float3 *dcos, float3 *normal, 
		float bulk_abscoeff, float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, float *Xc, float *Yc, struct histStruct * dPhotonHist,
		int boolToCollectData);
	#else
		int detection(float *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, 
		struct start_info info, float sensorRefl, float d_min, int* seed, float *dcos, float *normal, float bulk_abscoeff, float R, float xdetector, float ydetector, 
		unsigned long long int mynum_rebound, int ydim, int *h_num_detected_prim); 
	#endif

// calculate directional cosines of reflected/refracted vector.
	#ifdef USING_CUDA
		__device__ inline void trans_dir_cos(float3 *dcos, float3 *normal, float refl_theta, float trans_theta, int flag_ref, int mytid, struct start_info *info);  
	#else
		void trans_dir_cos(float *dcos, float *normal, float refl_theta, float trans_theta, int flag_ref, struct start_info info);
	#endif

// calculate new rough normal vector depending on value of 'beta' (roughness coefficient).
	#ifdef USING_CUDA
		__device__ inline void RoughSurface(float3 *normal, int2* seed, float beta);  
	#else
		void RoughSurface(float *normal, int* seed, float beta); 
	#endif



// RANECU pseudo random number generator
	#ifdef USING_CUDA
		__device__ inline void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int2* seed);	// initialize the generator
		__device__ inline int abMODm(int m, int a, int s);								// calculate a1*a2 MOD m
		__device__ inline float ranecu(int2* seed);									// Pseudo RNG returning float value (single-precision)
	#else
		void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int* seed);
		int abMODm(int m, int a, int s);
		float ranecu(int* seed);
	#endif


/////////////////////////////////////////
//
//       Global variables
//
/////////////////////////////////////////

	#ifdef USING_CUDA
		// counters
		__device__ unsigned long long int num_generated; // total # of photons generated for all the x-ray histories (across all threads)
		__device__ unsigned long long int num_detect;	 // total # of photons detected at the sensor plane of a column
		__device__ unsigned long long int num_abs_top;	 // total # of photons absorbed at the top surface of the detector (using 'top_absfrac')
		__device__ unsigned long long int num_abs_bulk;	 // total # of photons absorbed in the bulk of the detector (using 'bulk_abscoeff')
		__device__ unsigned long long int num_lost;	 // total # of photons lost when exiting out of the detector boundaries in x/y direction
		__device__ unsigned long long int num_outofcol;	 // total # of photons killed because they moved out of current column when reflected (due to precision errors)
		__device__ unsigned long long int num_theta1;	 // total # of photons killed if incidence angle > 1.57 or < 0 radian (after resampling max 100 times)	
		__device__ float photon_distance;     		 // total distance travelled by all the photons
		__device__ int dev_numPhotonHist;
	#else
		// counters
		unsigned long long int num_generated=0;	
		unsigned long long int num_detect=0;	
		unsigned long long int num_abs_top=0;	
		unsigned long long int num_abs_bulk=0;	
		unsigned long long int num_lost=0;	
		unsigned long long int num_outofcol=0;	 
		unsigned long long int local_counter=0;	 	// total number of photons terminated (either detected at sensor, absorbed at the top or in the bulk)
		unsigned long long int num_theta1=0;	

		//flags
		int absorbed=0;					// flag for photons absorbed at the top surface of the detector
		int detect=0;					// flag for photons detected at the sensor plane of the detector
		int bulk_abs=0;					// flag for photons absorbed in the bulk of a column

		float Xc=0.0f;					// center coordinates (x,y) of the current cylinder
		float Yc=0.0f;
		float photon_distance=0.0f;     		// total distance travelled by all the photons

		FILE *fp1;
	#endif

/////////////////////////////////////////
//
//    Functions definition
//
/////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// transports optical photon from its generation until it ends (detected/absorbed/lost).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
__global__ void algo(struct start_info *info, unsigned long long int *myimage, int *num_detected_primary, size_t pitch, int rowsread, float xdetector, float ydetector,  
float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x,  
float ubound_y, float d_max, float sensorRefl, struct histStruct * dPhotonHist, int boolToCollectData)
{
	unsigned long long int local_counter = 0; 	// total # optical photons terminated (either detected at sensor (bottom surface), absorbed at top or in the bulk). Used for checking if the total number of photons has been transported for this thread.
	unsigned long long int local_num_generated = 0;	// total # optical photons generated for this thread
	unsigned long long int num_rebound=0;
	float3 dcos, normal, pos; 			// directional cosine, normal and current position vector
	float rr=0.0f, theta=0.0f;
	float Xc=0.0f;					// center x,y coordinates of the current cylinder
	float Yc=0.0f;
	int tid = threadIdx.x + blockIdx.x * blockDim.x;	// cuda thread ID

	// flags
	int absorbed=0;			// flag for photons absorbed at the top surface of the detector (yes=1; no=0)
	int detect=0;			// flag for photons detected at the sensor plane of the detector
	int bulk_abs=0;			// flag for photons absorbed in the bulk of a column
	
	if(tid < rowsread)	// # threads launched is always > # scintillation events. Therefore, only the threads corresponding to an event works, rest do nothing. 
	{
		int NUM_EACH_THREAD = info[tid].str_N;		// number of photons to be simulated by each thread. info[tid] is the scintillation event for this thread.

		// Initialize vectors
		dcos.x = 0.0f; dcos.y = 0.0f; dcos.z = 0.0f;
		normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;
		pos.x = info[tid].str_x; 
		pos.y = info[tid].str_y; 
		pos.z = info[tid].str_z;	// starting location for this scintillation event.
		
		/*if(tid == 6 || tid == 7)
		{
				printf("algo - tid: %d, initial starting xyz %f %f %f\n", tid, pos.x, pos.y, pos.z);
		}*/
		
		int seed_input = 271828182+tid; 		// RANECU rng seed input
		int2 seed;

		// Initialize the RANECU generator in a position far away from the previous history
		init_PRNG(tid, 50000, seed_input, &seed);     

		// Initalize the device memory for directional cosine
		dcos.z = (ranecu(&seed) * 2.0f) - 1.0f;	// generate random number between -1.0 and 1.0

		rr = sqrt(1.0f - dcos.z*dcos.z);
		theta = ranecu(&seed) * twopipen;	// generate random number between 0 and 2pi
	
		dcos.x = rr*cos(theta);
		dcos.y = rr*sin(theta);

		// normalize
		if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
	 	{
			dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
			dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
		}

		local_num_generated++;		// increment number of generated photons

		/*if(tid == 0 || tid == 1 || tid == 3 || tid == 6 || tid == 7)
		{
			printf("tid: %d before while, local_num_generated: %lld, NUM_EACH_THREAD: %d\n", tid, local_num_generated, NUM_EACH_THREAD);
		}*/
				
		while(local_num_generated < (NUM_EACH_THREAD+1))	//run until NUM_EACH_THREAD particles generated
		{	
			/*if(tid == 6 || tid == 7)
			{
				printf("tid: %d after while(local_num_generated <...\n", tid);
			}*/
		
			if(absorbed == 0)  	  // not absorbed at the top surface of detector, check for absorption in the bulk and detection    
			{
				/*if(tid == 6 || tid == 7)
				{
					printf("tid: %d calling isotropic\n", tid);
				}*/
				
				bulk_abs = isotropic(&pos, &dcos, &seed, bulk_abscoeff, R, H, xdetector, ydetector, &info[tid], num_rebound, &Xc, &Yc, tid, dPhotonHist, boolToCollectData);
				if(bulk_abs == 0)	// not absorbed in the bulk, check for detection
				{
					/*if(tid == 6 || tid == 7)
					{
						printf("tid: %d calling detection\n", tid);
					}*/
				
					detect = detection(&pos, H, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, tid, num_detected_primary, 
					sensorRefl, d_min, &seed, &dcos, &normal, bulk_abscoeff, R, xdetector, ydetector, num_rebound, &Xc, &Yc, dPhotonHist, boolToCollectData);
				}
			}
		 
			if( ((detect == 1) || (absorbed == 1) || (bulk_abs == 1)) && (local_counter < (NUM_EACH_THREAD-1)) ) // particle terminated
			{
				/*if(tid == 6 || tid == 7)
				{
					printf("tid: %d particle terminated\n", tid);
				}*/
				
				local_counter++;	// increment # photons terminated

				// re-initialize all the arrays
				dcos.z = (ranecu(&seed) * 2.0f) - 1.0f;	

				rr = sqrt(1.0f - dcos.z*dcos.z);
				theta = ranecu(&seed) * twopipen;	
	
				dcos.x = rr*cos(theta);
				dcos.y = rr*sin(theta);

				if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
				{
					dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
					dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
					dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				}

				pos.x = info[tid].str_x; pos.y = info[tid].str_y; pos.z = info[tid].str_z; // set starting location of photon
				normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;
	
				if(beta > 0.0f)
					RoughSurface(&normal, &seed, beta);	// perturb smooth normal according to 'beta'

				local_num_generated++;
				absorbed = 0;
				detect = 0;
				bulk_abs = 0;
				num_rebound = 0;
				Xc = 0.0f;
				Yc = 0.0f;

			}
			else if( ((detect == 1) || (absorbed == 1) || (bulk_abs == 1)) && (local_counter == (NUM_EACH_THREAD-1)) )  // all the photons transported for this thread
			{
				/*if(tid == 6 || tid == 7)
				{
					printf("tid: %d all photons transported\n", tid);
				}*/
				local_counter++;
				break;
			}
			else if( (detect == 0) && (absorbed == 0) && (bulk_abs == 0) && (fabs(dcos.z - 0.0f) < epsilon) )  // check for trapped particle going back & forth dcos(z)=0
			{
				/*if(tid == 6 || tid == 7)
				{
					printf("tid: %d ccehck for trapped particles\n", tid);
				}*/
				
				// kill the particle and generate a new one instead - do not increment the local_counter
				// re-initialize all the arrays
		 		dcos.z = (ranecu(&seed) * 2.0f) - 1.0f;	

				rr = sqrt(1.0f - dcos.z*dcos.z);
				theta = ranecu(&seed) * twopipen;	
	
				dcos.x = rr*cos(theta);
				dcos.y = rr*sin(theta);

				if (((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) < (1.0f - epsilon)) || ((sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z)) > (1.0f + epsilon)))
				{
					dcos.x = dcos.x/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
					dcos.y = dcos.y/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
					dcos.z = dcos.z/(sqrt(dcos.x*dcos.x + dcos.y*dcos.y + dcos.z*dcos.z));
				}

				pos.x = info[tid].str_x; pos.y = info[tid].str_y; pos.z = info[tid].str_z;
				normal.x = 0.0f; normal.y = 0.0f; normal.z = 0.0f;
	
				if(beta > 0.0f)
					RoughSurface(&normal, &seed, beta);

				absorbed = 0;
				detect = 0;
				bulk_abs = 0;
				num_rebound = 0;
				Xc = 0.0f;
				Yc = 0.0f;
			}
			else	// transport photon
			{
				/*if(tid == 6 || tid == 7)
				{
					printf("tid: %d transport photon\n", tid);
				}*/
				num_rebound++;	// increment number of times the photon is rebounded from surface walls of a column
			    absorbed = boundary_analysis(&normal, &pos, &dcos, &seed, xdetector, ydetector, R, H, n1, n2, top_absfrac, beta, 
			    d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, &Xc, &Yc, pitch, info, tid, num_detected_primary, 
			    d_max, sensorRefl, dPhotonHist, boolToCollectData);
			}
		} 	// while loop ends

		atomicAdd(&num_generated, local_num_generated);

	}	// if tid=rowsread condition ends
	return;
}
#else		// C code

	int algo(float *normal, float *old_pos, float *pos, float *dcos, unsigned long long int *num_rebound, int* seed, struct start_info info, unsigned long long int *myimage, 
	float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, float bulk_abscoeff, float beta, float d_min, int pixelsize, float lbound_x, 
	float lbound_y, float ubound_x, float ubound_y, float sensorRefl, float d_max, int ydim, int *h_num_detected_prim)
	{

		float rr=0.0f, theta=0.0f, norm = 0.0f;		// used in calculating directional cosines
		float rnd_num = 0.0f;
		int myresult = 0;				// flag to check if the photon is terminated (yes=1; no=0)

		if(absorbed == 0)   		// not absorbed at the top surface of detector, check for absorption in the bulk and detection       
		 {
			bulk_abs = isotropic(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, info, num_rebound[local_counter]);

			if(bulk_abs == 0)	// not absorbed in the bulk, check for detection
			{
				detect = detection(pos, H, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, sensorRefl, d_min, seed, dcos, 
				normal, bulk_abscoeff, R, xdetector, ydetector, num_rebound[local_counter], ydim, h_num_detected_prim);
			}
		 }

		 
		if( (detect == 1) || (absorbed == 1) || (bulk_abs == 1) )	// photon terminated
		 {
			local_counter++;				// increment the photon counter

			// calculate directional cosines
			rnd_num = (ranecu(seed) * 2.0f) - 1.0f; 	// generate random number between (-1,1)	 	

			dcos[2] = rnd_num;				// generate random number between (-1,1)
			rr = sqrt(1.0-rnd_num*rnd_num);
			theta=ranecu(seed)*twopipen;			// generate random number between 0 and 2pi
			dcos[0]=rr*cos(theta);
			dcos[1]=rr*sin(theta);

			// normalize
			norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

			if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
			 {
				dcos[0] = dcos[0]/norm;
				dcos[1] = dcos[1]/norm;
				dcos[2] = dcos[2]/norm;
			 }


			pos[0] = info.str_x; pos[1] = info.str_y; pos[2] = info.str_z;	// set starting location of photon based on the PENELOPE scintillation events buffer
			old_pos[0] = info.str_x; old_pos[1] = info.str_y; old_pos[2] = info.str_z;
			normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 0.0f;		// initialize normal vector
	
			if(beta > 0.0f)
				RoughSurface(normal, seed, beta);	// perturb smooth normal according to 'beta' 

			absorbed = 0;
			detect = 0;
			bulk_abs = 0;

			myresult = 1;		// flag for photon termination

		 }
		else if( (detect == 0) && (absorbed == 0) && (bulk_abs == 0) && (fabs(dcos[2] - 0.0f) < epsilon) )  // check for trapped photon going back and forth dcos(z)=0
		 {
			// kill the photon and generate a new one instead - do not increment the counter

			// re-compute directional cosines
			rnd_num = (ranecu(seed) * 2.0f) - 1.0f; 
		 	
			dcos[2] = rnd_num;		
			rr = sqrt(1.0-rnd_num*rnd_num);
			theta=ranecu(seed)*twopipen;
			dcos[0]=rr*cos(theta);
			dcos[1]=rr*sin(theta);

			norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

			if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
			 {
				dcos[0] = dcos[0]/norm;
				dcos[1] = dcos[1]/norm;
				dcos[2] = dcos[2]/norm;
			 }

			pos[0] = info.str_x; pos[1] = info.str_y; pos[2] = info.str_z;
			normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 0.0f;
	
			if(beta > 0.0f)
				RoughSurface(normal, seed, beta);	

			absorbed = 0;
			detect = 0;
			bulk_abs = 0;

			myresult = 0;	// photon still alive
		 }
		else
		 {
			num_rebound[local_counter]++;	// increment the number of rebounds this photon undergoes from columnar walls
		    	absorbed = boundary_analysis(normal, pos, dcos, seed, xdetector, ydetector, R, H, n1, n2, top_absfrac, beta, d_min, 
		    	pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, h_num_detected_prim);

			myresult = 0;	// photon still alive
		 }

	return myresult;

	}

#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// determine where photon hits next within the column or if it gets absorbed in the material
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

	__device__ inline int isotropic(float3 *pos, float3 *dcos, int2* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector,
 	struct start_info *info, unsigned long long int mynum_rebound, float *Xc, float *Yc, int mytid, struct histStruct * dPhotonHist, int boolToCollectData)
	{
		float dsurf = 999.0f;		// distance to surface walls
		float dabs = 999.0f;		// distance to bulk absorption
		int flag_bulkabs = 0;		// flag to check if photon absorbed in the bulk
		
		dsurf = dist_to_surface(pos, dcos, R, H, xdetector, ydetector, info, mynum_rebound, Xc, Yc, mytid, dPhotonHist, boolToCollectData);	// distance to surface
		
		if (bulk_abscoeff > 0.0f)	
			dabs = (-1.0f/bulk_abscoeff) * log(ranecu(seed));					// distance to bulk absorption
		else
			dabs = 999999.0f;

		if (fabs(dsurf-(-99.0f)) < epsilon)				// photon lost: goes out of detector boundaries in dist_to_surface() function
		{
			flag_bulkabs = 1;
		}
		else if ( (dsurf < dabs) && (dsurf >= 0.0f) )			// photon not absorbed
		 {		
			flag_bulkabs = 0;
		 }
		else if ( (dsurf >= dabs) && (dabs >= 0.0f) )			// photon absorbed
		 {
			flag_bulkabs = 1;

			atomicAdd(&num_abs_bulk,1);				// increment to the number of optical photons absorbed in the bulk
			
			if(mytid < dev_numPhotonHist && boolToCollectData == 1)
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 4;
				dPhotonHist[mytid].histCounter ++;
			}
		 }

	   return flag_bulkabs;
	}

#else	// C code

	int isotropic(float *pos, float *dcos, int* seed, float bulk_abscoeff, float R, float H, float xdetector, float ydetector, 
	struct start_info info, unsigned long long int mynum_rebound)
	{
		float dsurf = 999.0f;		// distance to surface
		float dabs = 999.0f;		// distance to bulk absorption
		int flag_bulkabs = 0;		// flag to check if photon absorbed in the bulk

		dsurf = dist_to_surface(pos, dcos, R, H, xdetector, ydetector, info, mynum_rebound);	// distance to surface

		if (bulk_abscoeff > 0.0f)	
			dabs = (-1.0f/bulk_abscoeff) * log(ranecu(seed));				// distance to absorption
		else
			dabs = 999999.0f;

		if (fabs(dsurf-(-99.0f)) < epsilon)				// photon lost: goes out of detector boundaries in dist_to_surface() function
		{
			flag_bulkabs = 1;
		}
		else if ( (dsurf < dabs) && (dsurf >= 0.0f) )			// photon not absorbed
		 {		
			flag_bulkabs = 0;
		 }
		else if ( (dsurf >= dabs) && (dsurf >= 0.0f) )			// photon absorbed
		 {
			flag_bulkabs = 1;

			num_abs_bulk++;						// increment to the number of optical photons absorbed in the bulk
		 }

	   return flag_bulkabs;
	}
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate distance to surface (within same column)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

	__device__ float dist_to_surface(float3 *pos, float3 *dcos, float R, float H, float xdetector, float ydetector, struct start_info *info, 
	unsigned long long int mynum_rebound, float *Xc, float *Yc, int mytid, struct histStruct * dPhotonHist, int boolToCollectData)
	{
		float d=999.0f;				// distance to surface walls
		float d1=999.0f, d2=999.0f;
		float d_plane=999.0f, d_cyl=999.0f;	// distance to z-plane; to surface walls of the cylinder column
		float3 temp_pos = {0.0f};
		float3 my1 = {0.0f};
		float R1 = 999.0f, R2 = 999.0f;
		float stepsize = 0.5f;			// used for moving the photon towards the column, if it goes outside due to precision errors.
		int repeat = 0, ctr1 = 0; 		// number of times photon should be moved in steps towards the column before killing. In this code maximum steps = 100.
							// Valid only when goes out of column.
		
		/*if(mytid == 6 || mytid == 7)
		{
			printf("tid: %d in dist_to_surface\n", mytid);
		}*/
		
		// store current position in a temporary variable
		temp_pos.x = pos->x;
		temp_pos.y = pos->y;
		temp_pos.z = pos->z;

		// center of first column (assumed as x,y coordinates of the energy deposition event from Penelope)
		if(mynum_rebound == 0)				
		{
			*Xc = info->str_x;
			*Yc = info->str_y;
		}

		// solving quadratic equation for distance from a point to the surface of cylinder
		my1.x = dcos->x*dcos->x + dcos->y*dcos->y;
		my1.y = 2.0f*( ((pos->x-(*Xc))*dcos->x) + ((pos->y-(*Yc))*dcos->y) );
		my1.z = (pos->x-(*Xc))*(pos->x-(*Xc)) + (pos->y-(*Yc))*(pos->y-(*Yc)) - (R*R);
	
		// actual distance d = (d1 or d2)/sin_theta2	
		d1 = (-my1.y + (sqrt( (my1.y*my1.y) - (4.0f*my1.x*my1.z) )))/(2.0f * my1.x);
		d2 = (-my1.y - (sqrt( (my1.y*my1.y) - (4.0f*my1.x*my1.z) )))/(2.0f * my1.x);


		// hits either upper half surface or top of the cylinder
		if(dcos->z > 0.0f) 
		 {

		  if((fabs(dcos->x - 0.0f) < epsilon) && (fabs(dcos->y - 0.0f) < epsilon) && (fabs(dcos->z - 1.0f) < epsilon))  
		  // if photon travel straight in +z axis direction
		   {
			d = (H/2.0f - pos->z)/dcos->z;

			pos->z = H/2.0f;	

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;
		   }
		  else
		   {		

			// calculate distance to infinite plane at z=H/2
			d_plane = (H/2.0f - pos->z)/(dcos->z);

			// calculate the distance to the upper half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos->z = temp_pos.z + d*dcos->z;
			 }
			else
			 {
				d = d_plane;		
				pos->z = H/2.0f;
			 }

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		   }	// else loop ends
	
		
		 } // if loop for dcos.z > 0 ends

		else if(dcos->z < 0.0f) // hits either lower half or bottom of cylinder
		 {
		  // if photon travels in -Z direction staright, then it should get detected
		  if ((fabs(dcos->x-0.0f) < epsilon) && (fabs(dcos->y-0.0f) < epsilon) && (fabs(dcos->z - (-1.0f)) < epsilon))  
		   {
			d = (-H/2.0f - pos->z)/dcos->z;  
			pos->z = -H/2.0f;
				
			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;
		   }
		  else
		   {

			// calculate distance to infinite plane at z=-H/2
			d_plane = (-H/2.0f - pos->z)/(dcos->z);

			// calculate the distance to the lower half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos->z = temp_pos.z + d*dcos->z;
			 }
			else
			 {
				d = d_plane;		
				pos->z = -H/2.0f;
			 }

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		   }	// else loop ends

		 }	// else if loop for dcos.z < 0 ends

		else	// when dcos.z=0.0 (will hit only the side of the cylinder)
		 {
			// calculate the distance to the side of cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		
		
			d = d_cyl;
			pos->z = temp_pos.z + d*dcos->z;

			pos->x = temp_pos.x + d*dcos->x;
			pos->y = temp_pos.y + d*dcos->y;

		 }

		// condition to check that pos is within detector boundaries - if true, photon LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) || (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
			{
				d = -99.0f;
				atomicAdd(&num_lost,1);
				goto distexit;
			}
		else
			atomicAdd(&photon_distance, d);		// add distance travelled to global variable


		// check if photon is outside the current column. This can happen due to single precision errors. If yes, then move the photon towards the column based on the 'stepsize' and check again. Repeat this 100 times, if still outside then kill the photon.
		R1 = sqrt((pos->x - (*Xc))*(pos->x - (*Xc)) + (pos->y - (*Yc))*(pos->y - (*Yc)));

		repeat = 0;	// counters for checking how many times has the photon been moved towards the column.
		ctr1 = 0;

		while( (R1 > (R-1e-5)) && (repeat < 10) && (ctr1 < 10) ) // R1 > R-some small value; to avoid single precision errors
		{

			// store current position
			temp_pos.x = pos->x;
			temp_pos.y = pos->y;

			// move photon by stepsize (in microns) in the incident direction
			pos->x = pos->x + stepsize*(-dcos->x);
			pos->y = pos->y + stepsize*(-dcos->y);

			R2 = sqrt((pos->x - (*Xc))*(pos->x - (*Xc)) + (pos->y - (*Yc))*(pos->y - (*Yc))); // new radius

			if(R2 > R1) // means the photon is moving farther away from the column 
			{	    // this can happen if the stepsize is big enough for the photon to pass through the column and exit on the other side.

				// move it back to previous position, reduce the stepsize and try moving it again.
				pos->x = temp_pos.x;
				pos->y = temp_pos.y;

				stepsize = stepsize/2.0f;
				ctr1++;
			}
			else
			{
				R1 = R2;
				repeat++;
			}

		}

		// kill the photon if still outside the column
		if(R1 > (R-1e-5))
		 {
			d = -99.0f;
			atomicAdd(&num_outofcol,1);
			goto distexit;
		 }

		// condition to check that pos is within detector boundaries - if true, photon LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) || (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
			{
				d = -99.0f;
				atomicAdd(&num_lost,1);
				goto distexit;
			}

	distexit:
		if(mytid < dev_numPhotonHist && boolToCollectData == 1)
		{
			if(d == -99.0f)
			{
				//photon lost
			}
			else
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
				dPhotonHist[mytid].histCounter ++;
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("dist_to_surface - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
			}
		}
		return d;
	}	

#else	// C code

	float dist_to_surface(float *pos, float *dcos, float R, float H, float xdetector, float ydetector, struct start_info info, 
	unsigned long long int mynum_rebound)
	{
		float d=999.0f;				// distance to surface walls
		float d1=999.0f, d2=999.0f;
		float d_plane=999.0f, d_cyl=999.0f;	// distance to z-plane; to surface walls of the cylinder column
		float temp_pos[3] = {0.0f};
		float my1[3] = {0.0f};
		float R1 = 999.0f, R2 = 999.0f;
		float stepsize = 0.5f;			// used for moving the photon towards the column, if it goes outside due to precision errors.
		int repeat = 0, ctr1 = 0; 		// number of times photon should be moved in steps towards the column before killing. In this code maximum steps = 100.
							// Valid only when goes out of column.

		// store current position in a temporary variable
		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		// center of first column (assumed as x,y coordinates of the energy deposition event from Penelope)
		if(mynum_rebound == 0)				
		{
			Xc = info.str_x;
			Yc = info.str_y;
		}

		// solving quadratic equation for distance from a point to the surface of cylinder
		my1[0] = dcos[0]*dcos[0] + dcos[1]*dcos[1];
		my1[1] = 2.0f*( ((pos[0]-(Xc))*dcos[0]) + ((pos[1]-(Yc))*dcos[1]) );
		my1[2] = (pos[0]-(Xc))*(pos[0]-(Xc)) + (pos[1]-(Yc))*(pos[1]-(Yc)) - (R*R);
	
		// actual distance d = (d1 or d2)/sin_theta2	
		d1 = (-my1[1] + (sqrt( (my1[1]*my1[1]) - (4.0f*my1[0]*my1[2]) )))/(2.0f * my1[0]);
		d2 = (-my1[1] - (sqrt( (my1[1]*my1[1]) - (4.0f*my1[0]*my1[2]) )))/(2.0f * my1[0]);

	
		// hits either upper half surface or top surface of cylinder
		if(dcos[2] > 0.0f) 
		 {

		  if((fabs(dcos[0] - 0.0f) < epsilon) && (fabs(dcos[1] - 0.0f) < epsilon) && (fabs(dcos[2] - 1.0f) < epsilon))  
		  // if photon travel straight in +z axis direction
		   {
			d = (H/2.0f - pos[2])/dcos[2];

			pos[2] = H/2.0f;	

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];
		   }
		  else
		   {		

			// calculate distance to infinite plane at z=H/2
			d_plane = (H/2.0f - pos[2])/(dcos[2]);

			// calculate the distance to the upper half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos[2] = temp_pos[2] + d*dcos[2];
			 }
			else
			 {
				d = d_plane;		
				pos[2] = H/2.0f;
			 }

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		   }	// else loop ends
	
		
		 } // if loop for dcos[2] > 0 ends

		else if(dcos[2] < 0.0f) // hits either lower half or bottom of cylinder
		 {
		  // if photon travels in -Z direction staright, then it should get detected
		  if ((fabs(dcos[0]-0.0f) < epsilon) && (fabs(dcos[1]-0.0f) < epsilon) && (fabs(dcos[2] - (-1.0f)) < epsilon))  
		   {
			d = (-H/2.0f - pos[2])/dcos[2];  
			pos[2] = -H/2.0f;
				
			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];
		   }
		  else
		   {

			// calculate distance to infinite plane at z=-H/2
			d_plane = (-H/2.0f - pos[2])/(dcos[2]);

			// calculate the distance to the lower half of the cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		

			// find min from d_plane and d_cyl
			if(d_plane >= d_cyl)
			 {
				d = d_cyl;
				pos[2] = temp_pos[2] + d*dcos[2];
			 }
			else
			 {
				d = d_plane;		
				pos[2] = -H/2.0f;
			 }

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		   }	// else loop ends

		 }	// else if loop for dcos[2] < 0 ends

		else	// when dcos[2]=0.0 (will hit only the side of the cylinder)
		 {
			// calculate the distance to the side of cylinder
			if(d1 >= d2) 
			{
				d_cyl = d1;
			}
			else if(d2 > d1)
			{
				d_cyl = d2;
			}
		
		
			d = d_cyl;
			pos[2] = temp_pos[2] + d*dcos[2];

			pos[0] = temp_pos[0] + d*dcos[0];
			pos[1] = temp_pos[1] + d*dcos[1];

		 }

		// condition to check that pos is within detector boundaries - if true, photon LOST
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
			{
				d = -99.0f;
				num_lost++;
				goto distexit;
			}
		else
			photon_distance = photon_distance + d;		// add distance travelled to global variable


		// check if photon is outside of current column. This can happen due to single precision errors. If yes, then move the photon towards the column based on the 'stepsize' and check again. Repeat this 100 times, if still outside then kill the photon.
		R1 = sqrt((pos[0] - Xc)*(pos[0] - Xc) + (pos[1] - Yc)*(pos[1] - Yc));
		
		repeat = 0;
		ctr1 = 0;

		while( (R1 > (R-1e-5)) && (repeat < 10) && (ctr1 < 10) ) // R1 > R1-some small value..because of single precision errors that comparison with R may generate
		{

			// store current position
			temp_pos[0] = pos[0];
			temp_pos[1] = pos[1];

			// move photon by 0.5 um in the incident direction
			pos[0] = pos[0] + stepsize*(-dcos[0]);
			pos[1] = pos[1] + stepsize*(-dcos[1]);

			R2 = sqrt((pos[0] - Xc)*(pos[0] - Xc) + (pos[1] - Yc)*(pos[1] - Yc));

			if(R2 > R1) // means the photon is moving farther away from the column 
			{	    // this can happen if the stepsize is big enough for the photon to pass through the column and exit on other side.

				// move it back to previous position, reduce the stepsize and try moving it again.
				pos[0] = temp_pos[0];
				pos[1] = temp_pos[1];

				stepsize = stepsize/2.0f;
				ctr1++;
			}
			else
			{
				R1 = R2;
				repeat++;
			}

		}

		// kill the photon if still outside the column
		if(R1 > (R-1e-5))
		 {
			d = -99.0f;
			num_outofcol++;
			goto distexit;
		 }

		// condition to check that pos is within detector boundaries - if true, photon LOST
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) || (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
			{
				d = -99.0f;
				num_lost++;
				goto distexit;
			}

	distexit:
	 return d;
	}
	
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate the directional cosines of the reflected vector
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
__device__ int boundary_analysis(float3 *normal, float3 *pos, float3 *dcos, int2* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, 
float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, 
float *Xc, float *Yc, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float d_max, float sensorRefl, struct histStruct * dPhotonHist,
int boolToCollectData)
{
	float3 dcos_temp = {0.0f};
	float Pr = 0.0f, Pt = 0.0f;	// Probability of reflection and transmission
	float theta1 = 0.0f;		// angle between normal and reflected (radians)
	float theta2 = 0.0f;		// angle between normal and transmitted (radians)
	float cct1 = 0.0f;		// columnar crosstalk
	int trans_flag = 0.0f;		// flag - photon terminated during transmission to a new column; used in transmit() (terminated: yes=1; no=0)
	int flag_abs = 0;		// flag - photon absorbed at the top surface (yes=1; no=0)
	int flag_call_transmit = 1;	// flag - photon moves within a column (flag = 0) [call isotropic()] or between columns (flag = 1) [call transmit()]
	int flagCCT = 0;		// flag - photon cross over (yes=1; no=0); send as input to transmit()
	int newnormalctr=0;		// counter - if angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate normal; max. 100 times)
	int theta1ctr=0;		// counter - if theta1 > 1.57 or < 0 radians (recalculate normal; max. 100 times)
	int oldN_Rctr=0;		// counter - if angle between reflected dir. cosine and smooth normal > 1.57 radians (recalculate max. 25 times)
	int reperturb_ctr = 0;		// counter - if angle between reflected dir. cosine and rough normal > 1.57 radians (reperturb normal; max. 3 times)
	float newdepth = 0.0f;		// bottom depth for which CCT=1 (z_a)
	float temp_norm = 0.0f;
	float mag = 0.0f;
	float rr_rnd = 0.0f, theta_rnd = 0.0f;
	float3 old_normal = {0.0f};
	float3 old_dcos = {0.0f};
	float angle_oldN_R = 0.0f;
	float cos_newangle = 0.0f, newangle = 0.0f;

	/*if(mytid == 6 || mytid == 7)
	{
			printf("tid: %d in boundary_analysis\n", mytid);
	}*/
	
	// determine the coordinates of normal
	if ( (fabs(pos->z - (float)(H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	// reached top surface and dir. cosine in z-direction is positive
	{
		if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )	// photon gets absorbed; 'top_absfrac' is the top surface absorption fraction (0,1)
		{
			flag_abs = 1;
			atomicAdd(&num_abs_top, 1);				// increment # photons absorbed at the top

			if(mytid < dev_numPhotonHist && boolToCollectData == 1)
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 3;
				dPhotonHist[mytid].histCounter ++;
			}
			
			/*if(mytid == 6 || mytid == 7)
			{
				printf("boundary analysis - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
			}*/
		}
		else
		{
			// normal		
			normal->x = 0.0f;
			normal->y = 0.0f;
			normal->z = -1.0f;

			// assign new directional cosines - top surface being isotropic reflector
			dcos->z = -fabs((ranecu(seed) * 2.0f) - 1.0f);
			rr_rnd = sqrt(1.0f - dcos->z*dcos->z);
			theta_rnd = ranecu(seed)*twopipen;	

			dcos->x=rr_rnd*cos(theta_rnd);
			dcos->y=rr_rnd*sin(theta_rnd);
	
			flag_abs = 0;
		}
	}	
	else 	// gets reflected or transmitted
	{	
		// Columnar crosstalk
		newdepth = H*0.2f;	// top 20% depth CCT=1. considering CsI layer only. NO organic polymer coating.
			
		if( (pos->z <= H/2.0f) && (pos->z >= (H/2.0f - newdepth)) )	// top 20% - 100% cct
		{
			cct1 = 1.0f;
		}
		else if( (pos->z < (H/2.0f - newdepth)) && (pos->z >= 0.0f) )  // from 20% depth to 50% - linear 100% to 50% 
		{
			cct1 = (pos->z/(2.0f*(H/2.0f - newdepth))) + 0.5;	
		}
		else if( (pos->z < 0.0f) && (pos->z >= (-H/2.0f)) ) // bottom 50% to -H/2 - 50% to 100% CCT
		{
			cct1 = ( (pos->z - (-H/2.0f))/(2.0f * (-H/2.0f)) ) + 1.0 ;
		}
	
		if(ranecu(seed) < cct1)		// columnar crosstalk occurs
		{
			// photon crosses over to adjacent column with random orientation. directional cosine do not change.
			flagCCT = 1;

			trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, 
			lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, 
			flagCCT, Xc, Yc, dPhotonHist, boolToCollectData);

			if (trans_flag == 1)		// photon terminated during transmission
				flag_abs = 1;
			else if (trans_flag == 0)	// photon still alive; transmitted to adjacent column
			{
				// calculate new column's center coordinates
				*Xc = (float)( pos->x + R*(-normal->x) );
				*Yc = (float)( pos->y + R*(-normal->y) );

				if(mytid < dev_numPhotonHist && boolToCollectData == 1)
				{
					dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
					dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
					dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
					dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
					dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
					dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
					dPhotonHist[mytid].histCounter ++;
				}
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("boundary analysis - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
				flag_abs = 0;
			}
		}
		else	// no cross over
		{
			prpt:
				// within the column
				if(flag_call_transmit == 1)			// photon is currrently within a column with center Xc,Yc
				{
					mag = sqrt( (((*Xc)-pos->x) * ((*Xc)-pos->x)) + (((*Yc)-pos->y) * ((*Yc)-pos->y)) );
					normal->x = ((*Xc)-pos->x)/mag;
					normal->y = ((*Yc)-pos->y)/mag;
					normal->z = 0.0f;
		
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	// perturb normal for rough surface

					flag_abs = 0;
				}
				// outside the column
				else if (flag_call_transmit == 0)		// photon is currently between columns and has not entered any column yet. New normal is sampled in the transmit(), so do not calculate normal here.
				{
					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					*Xc = (float)( pos->x + R*(-normal->x) );
					*Yc = (float)( pos->y + R*(-normal->y) );

					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("boundary analysis - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
			
				    flag_abs = 0;
				}

				
				dcos_temp.x = -dcos->x;		// -dcos -> invert the incident dcos vector; to get the smaller angle between normal and dcos
				dcos_temp.y = -dcos->y;
				dcos_temp.z = -dcos->z;

				old_normal.x = normal->x;
				old_normal.y = normal->y;
				old_normal.z = normal->z;
	
				old_dcos.x = dcos->x;
				old_dcos.y = dcos->y;
				old_dcos.z = dcos->z;

			reperturb:
				normal->x = old_normal.x;
				normal->y = old_normal.y;
				normal->z = old_normal.z;

				dcos->x = old_dcos.x;
				dcos->y = old_dcos.y;
				dcos->z = old_dcos.z;

				dcos_temp.x = -dcos->x;
				dcos_temp.y = -dcos->y;
				dcos_temp.z = -dcos->z;

				if( (flag_call_transmit == 1) && (reperturb_ctr != 0) )		// within the column
				{
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	
				}
				if( (flag_call_transmit == 0) && (reperturb_ctr != 0) )		// outside the column
				{
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	

					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					*Xc = (float)( pos->x + R*(-normal->x) );
					*Yc = (float)( pos->y + R*(-normal->y) );

					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("boundary analysis - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
				}

				// Using Snell's law, calculate theta1 (angle between normal and reflected) and theta2 (angle between normal and transmitted)
			no_perturbation:
				theta1 = dot_product(&dcos_temp,normal);		// cosine of angle between incident in opposite direction and normal (in radians)

				if ( (theta1 > 1.0f) || (theta1 < 0.0f) )	// if incidence angle > 1.57 radian or < 0 radian, then recalculate normal
				{
					// if photon was transmitted, then new normal has to be sampled again
					if(flag_call_transmit == 0)
					{
					mynewnormal:
						normal->x = dcos_temp.x;		// normal = inverted dcos of incident vector
						normal->y = dcos_temp.y;
						normal->z = dcos_temp.z;

						RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

						mag = sqrt(normal->x*normal->x + normal->y*normal->y);

						// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
						normal->z = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
						normal->x = normal->x/mag;		
						normal->y = normal->y/mag;

						if(beta > 0.0f)
							RoughSurface(normal, seed, beta);

						// find the angle between normal and -dcos
						cos_newangle = dot_product(&dcos_temp, normal);
						newangle = acosf(cos_newangle);

						if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new rough normal is within +- 1.57 radians from inverted dcos
						{						// keep looping until get a theta within 1.57 radians - maximum iterations 100
							if(newnormalctr < 100)
							{
								newnormalctr++;
								goto mynewnormal;			
							}
							else 					// else terminate the photon
							{
								atomicAdd(&num_theta1,1);  // increment the counter - # photons terminated due to incidence angle > 1.57 or < 0 radian
								flag_abs = 1;
								newnormalctr = 0;
								goto baexit;
							}
						}

					}
		
					if(theta1ctr < 100)	// recalculate max. 100 times	
					{
						theta1ctr++;
						goto prpt;
					}
					else			// terminate photon
					{
						atomicAdd(&num_theta1,1);	
						flag_abs = 1;
						theta1ctr = 0;
						goto baexit;
					}
				}
				else				// 0 < theta1 < 1.57 (radian); continue the photon transport
					theta1 = acosf(theta1);
		

				// check for conditions where photon can only reflect
				if (flag_call_transmit == 1)	// only valid when photon is within the column and can transmit outside the column. asin(n1/n2) -> nan
				{
					if (theta1 > asin(n2/n1))	// critical angle condition for total internal reflection (TIR)
					{
						Pr = 1.0f;		// TIR occurs
						Pt = 0.0f;
					}
			       		else if ( theta1 < epsilon ) 	// theta1 ~= 0, always reflect
					{
				        	theta1 = 0.00042;       			// make theta1 = very small number, to avoid getting nan probabilities
					        theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians

					        // Using Fresnel's law, compute probability of reflection and transmission 
					        Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				        	Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else    // the ray will transmit
					{
				        	theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians
					        
					        // Using Fresnel's law, compute probability of reflection and transmission 
				        	Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
					        Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}
				else if (flag_call_transmit == 0)	// outside the column
				{		
					if((n1/n2) < 1.57f)	
					{
						if (theta1 > asin(n1/n2))	// critical angle condition for total internal reflection (TIR)
						{
							Pr = 1.0f;		// TIR occurs
							Pt = 0.0f;
						}
					}
					else if ( theta1 < epsilon )	// theta1 ~= 0, then always reflect
					{
						theta1 = 0.00042;	// make theta1 a very smal number, to avoid getting nan probabilities
						theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
	
						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else	// photon transmits
					{
						theta2 = asinf((float)(n2/n1)*sin(theta1));

						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
				}


				// normalize Pr and Pt
				temp_norm = Pr + Pt;
				Pr = Pr/temp_norm;
				Pt = Pt/temp_norm;

				if(ranecu(seed) < Pr)				// reflection
				{
					trans_dir_cos(dcos, normal, theta1, theta2, 0, mytid, info);	// compute new directional cosines

					// check that reflected vector is within 1.57 radians from original normal
					angle_oldN_R = dot_product(&old_normal, dcos);
					angle_oldN_R = acosf(angle_oldN_R);


					if (angle_oldN_R > 1.57f) 		// > 1.57 radians, reperturb the normal
					{
						reperturb_ctr++;

						if(reperturb_ctr < 4)		// reperturb maximum 3 times
							goto reperturb;
						else				// else calculate using smooth surface normal (old_normal)
						{
							normal->x = old_normal.x;
							normal->y = old_normal.y;
							normal->z = old_normal.z;

							dcos->x = old_dcos.x;
							dcos->y = old_dcos.y;
							dcos->z = old_dcos.z;

							dcos_temp.x = -dcos->x;
							dcos_temp.y = -dcos->y;
							dcos_temp.z = -dcos->z;

							reperturb_ctr = 0;
				
							if(oldN_Rctr < 25)	// resample max 75 times (25 * reperturb 3 times)	
							{
								oldN_Rctr++;
								goto no_perturbation;
							}
							else			// terminate photon
							{
								atomicAdd(&num_theta1,1);// increment the counter - # photons terminated due to incidence angle > 1.57 or < 0 radian
								flag_abs = 1;
								oldN_Rctr = 0;
								goto baexit;
							}
						}
					}

					if (flag_call_transmit == 0)		// reflects between columns; calculate distance using transmit()
					{
						flagCCT = 0;	

						trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, 
						lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, 
						flagCCT, Xc, Yc, dPhotonHist, boolToCollectData);

						if (trans_flag == 1)		// photon terminated
							flag_abs = 1;
						else if (trans_flag == 0)
							goto prpt;				
					}
				}
				else						// transmits 
				{
					trans_dir_cos(dcos, normal, theta1, theta2, 1, mytid, info);		// compute transmitted directional cosines

					if (flag_call_transmit == 1)		// photon travels between columns
					{
						flag_call_transmit = 0;
						flagCCT = 0;

						trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, 
						lbound_x, lbound_y, ubound_x, ubound_y, myimage, pitch, info, mytid, num_detected_primary, d_max, sensorRefl, 
						flagCCT, Xc, Yc, dPhotonHist, boolToCollectData);

						if (trans_flag == 1)		// photon terminated
							flag_abs = 1;
						else if (trans_flag == 0)	// hits a column
							goto prpt;		// check again to see if it gets reflected or transmitted
					}			
				}
		
			} // else 'prpt ends
		
		} // main else ends

	baexit:
	   return flag_abs;

	}	

#else	// C code

	int boundary_analysis(float *normal, float *pos, float *dcos, int* seed, float xdetector, float ydetector, float R, float H, float n1, float n2, float top_absfrac, float beta, 	float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, 
	float sensorRefl, int ydim, int *h_num_detected_prim)
	{

		float dcos_temp[3] = {0.0f};
		float Pr = 0.0f, Pt = 0.0f;	// Probability of reflection and transmission
		float theta1 = 0.0f;		// angle between normal and reflected vector (radians)
		float theta2 = 0.0f;		// angle between normal and transmitted vector (radians)
		float cct1 = 0.0f;		// columnar crosstalk
		int trans_flag = 0.0f;		// flag - photon terminated during transmission to a new column; used in transmit() (terminated: yes=1; no=0)
		int flag_abs = 0;		// flag - photon absorbed at the top surface (yes=1; no=0)
		int flag_call_transmit = 1;	// flag - photon moves within a column (flag = 0) [call isotropic()] or between columns (flag = 1) [call transmit()]
		int flagCCT = 0;		// flag - photon cross over (yes=1; no=0); send as input to transmit()
		int newnormalctr=0;		// counter - if angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate normal; max. 100 times)
		int theta1ctr=0;		// counter - if theta1 > 1.57 or < 0 radians (recalculate normal; max. 100 times)
		int oldN_Rctr=0;		// counter - if angle between reflected dir. cosine and smooth normal > 1.57 radians (recalculate max. 25 times)
		int reperturb_ctr = 0;		// counter - if angle between reflected dir. cosine and rough normal > 1.57 radians (reperturb normal; max. 3 times)
		float newdepth = 0.0f;		// bottom depth for which CCT=1 (z_a)
		float temp_norm = 0.0f;
		float mag = 0.0f;
		float rr_rnd = 0.0f, theta_rnd = 0.0f;
		float old_normal[3] = {0.0f};
		float old_dcos[3] = {0.0f};
		float angle_oldN_R = 0.0f;
		float cos_newangle = 0.0f, newangle = 0.0f;


		// determine the coordinates of normal
		if ( (fabs(pos[2] - (float)(H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	// reached top surface and dir. cosine in z-direction is positive
		{
	
			if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )	// photon gets absorbed; 'top_absfrac' is the top surface absorption fraction (0,1)
			{
				flag_abs = 1;
				num_abs_top++;						// increment # photons absorbed at the top surface
			}
			else
			{
				normal[0] = 0.0f;
				normal[1] = 0.0f;
				normal[2] = -1.0f;

				// assign new directional cosines; top surface is isotropic reflector
				dcos[2] = -fabs((ranecu(seed) * 2.0f) - 1.0f);
				rr_rnd = sqrt(1.0f - dcos[2]*dcos[2]);
				theta_rnd = ranecu(seed)*twopipen;	

				dcos[0]=rr_rnd*cos(theta_rnd);
				dcos[1]=rr_rnd*sin(theta_rnd);
	
				flag_abs = 0;
			}

		}	
		else 	// photon reflected or transmitted
		{	
			// Columnar crosstalk
			newdepth = H*0.2f;	// top 20% depth CCT=1. considering CsI layer only. NO organic polymer coating.
		
			if( (pos[2] <= H/2.0f) && (pos[2] >= (H/2.0f - newdepth)) )	// top 20% - 100% cct
			{
				cct1 = 1.0f;
			}
			else if( (pos[2] < (H/2.0f - newdepth)) && (pos[2] >= 0.0f) )  // from 20% depth to 50% - linear 100% to 50% 
			{
				cct1 = (pos[2]/(2.0f*(H/2.0f - newdepth))) + 0.5;	
			}
			else if( (pos[2] < 0.0f) && (pos[2] >= (-H/2.0f)) ) // bottom 50% to (-H/2 - 4 um polymer) - 50% to 100% CCT
			{
				cct1 = ( (pos[2] - (-H/2.0f))/(2.0f * (-H/2.0f)) ) + 1.0 ;
			}

	
			if(ranecu(seed) < cct1)		// columnar crosstalk occurs
			{

				// photon crosses over to adjacent column with random orientation. directional cosine do not change.
				flagCCT = 1;

				trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

				if (trans_flag == 1)		// photon terminated
					flag_abs = 1;
				else if (trans_flag == 0)	// photon still alive; crossed over to adjacent column
				{
					// calculate new column's center coordinates
					Xc = (float)( pos[0] + R*(-normal[0]) );
					Yc = (float)( pos[1] + R*(-normal[1]) );

					flag_abs = 0;
				}
			}
			else
			{

			prpt:

				// within the column
				if(flag_call_transmit == 1)			// photon is currently within a column with center Xc,Yc
				{
					mag = sqrt( (((Xc)-pos[0]) * ((Xc)-pos[0])) + (((Yc)-pos[1]) * ((Yc)-pos[1])) );
					normal[0] = ((Xc)-pos[0])/mag;
					normal[1] = ((Yc)-pos[1])/mag;
					normal[2] = 0.0f;
		
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	// perturb normal for rough surface according to 'beta'

					flag_abs = 0;
				}
				// outside the column
				else if (flag_call_transmit == 0)		// photon is currently between columns and has not entered any column yet. New normal is sampled in the transmit(), so do not calculate normal here.
				{
					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					Xc = (float)( pos[0] + R*(-normal[0]) );
					Yc = (float)( pos[1] + R*(-normal[1]) );

				        flag_abs = 0;
				}

				dcos_temp[0] = -dcos[0];	// -dcos -> invert the incident dcos vector; to get the smaller angle between normal and dcos
				dcos_temp[1] = -dcos[1];
				dcos_temp[2] = -dcos[2];

				old_normal[0] = normal[0];
				old_normal[1] = normal[1];
				old_normal[2] = normal[2];
	
				old_dcos[0] = dcos[0];
				old_dcos[1] = dcos[1];
				old_dcos[2] = dcos[2];

			reperturb:
				normal[0] = old_normal[0];
				normal[1] = old_normal[1];
				normal[2] = old_normal[2];

				dcos[0] = old_dcos[0];
				dcos[1] = old_dcos[1];
				dcos[2] = old_dcos[2];

				dcos_temp[0] = -dcos[0];
				dcos_temp[1] = -dcos[1];
				dcos_temp[2] = -dcos[2];

				if( (flag_call_transmit == 1) && (reperturb_ctr != 0) )		// within the column
				 {
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	
				 }
				if( (flag_call_transmit == 0) && (reperturb_ctr != 0) )		// outside the column
				 {
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);	

					// center of new column (obtained by inverting the new normal sampled in transmit() and finding center at distance R from current position)
					Xc = (float)( pos[0] + R*(-normal[0]) );
					Yc = (float)( pos[1] + R*(-normal[1]) );

				 }
			
				// Using Snell's law, calculate theta1 (angle between normal and reflected) and theta2 (angle between normal and transmitted)
			no_perturbation:
				theta1 = dot_product(dcos_temp,normal);		// cosine of angle between incident in opposite direction and normal (in radians)

		
				if ( (theta1 > 1.0f) || (theta1 < 0.0f) )	// if incidence angle > 1.57 or < 0 radians, then recalculate normal
				{
					// if photon was transmitted, then new normal has to be sampled again
					if(flag_call_transmit == 0)
					{
					mynewnormal:
						normal[0] = dcos_temp[0];		// normal = inverted dir. cosine of the incident vector
						normal[1] = dcos_temp[1];
						normal[2] = dcos_temp[2];

						RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

						mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

						// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
						normal[2] = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
						normal[0] = normal[0]/mag;		
						normal[1] = normal[1]/mag;

						if(beta > 0.0f)
							RoughSurface(normal, seed, beta);

						// find the angle between normal and -dcos
						cos_newangle = dot_product(dcos_temp, normal);
						newangle = acosf(cos_newangle);

						if ( (newangle < 0.0f) || (newangle > 1.57f) )	// new normal within +- 1.57 radians from inverted dcos
						{						// keep looping until get a theta within 1.57 radians - maximum iterations 100
							if(newnormalctr < 100)
							{
								newnormalctr++;
								goto mynewnormal;			
							}
							else 					// else terminate photon
							{
								num_theta1++;
								flag_abs = 1;
								newnormalctr = 0;
								goto baexit;
							}				
						}

					}


					if(theta1ctr < 100)	// recalculate max. 100 times
					{
						theta1ctr++;
						goto prpt;
					}
					else			// terminate photon
					{
						num_theta1++;	// increment the counter for 'theta1' - # photons terminated due to incidence angle > 1.57 or < 0 radian
						flag_abs = 1;
						theta1ctr = 0;
						goto baexit;
					}
				}
				else		// 0< theta1 < 1.57 radians; continue photon transport
					theta1 = acosf(theta1);
		

				// check for conditions where photon can only reflect
				if (flag_call_transmit == 1)	// only valid when photon within the column and can transmit outside the column. asin(n1/n2) -> nan
				{
					if (theta1 > asin(n2/n1))	// critical angle condition for total internal reflection (TIR)
					{
						Pr = 1.0f;		// TIR occurs
						Pt = 0.0f;
					}
			       		else if ( theta1 < epsilon ) 	// theta1 ~= 0, then always reflect
					{
				        	theta1 = 0.00042;       // make theta1 a very small number, to avoid getting nan probabilities
					        theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians

					        // Using Fresnel's law, compute probability of reflection and transmission 
					        Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
				        	Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else    // the ray will transmit
					{
				        	theta2 = asinf((float)(n1/n2)*sin(theta1));     // refracted/transmitted angle in radians
					        
					        // Using Fresnel's law, compute probability of reflection and transmission 
				        	Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
					        Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}
				else if (flag_call_transmit == 0)	// outside the column
				{		
					if((n1/n2) < 1.57f)	
					{
						if (theta1 > asin(n1/n2))	// critical angle condition for total internal reflection (TIR)
						{
							Pr = 1.0f;		// TIR occurs
							Pt = 0.0f;
						}
					}
					else if ( theta1 < epsilon )	// theta1 ~= 0, then always reflect
					{
						theta1 = 0.00042;	// make theta1 a very small number, to avoid getting nan probabilities
						theta2 = asinf((float)(n1/n2)*sin(theta1)); 	// refracted/transmitted angle in radians
	
						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}
					else	// photon transmits
					{
						theta2 = asinf((float)(n2/n1)*sin(theta1));

						// Using Fresnel's law, compute probability of reflection and transmission 
						Pr = 1/2.0f*( (pow(tan(theta1-theta2),2)/pow(tan(theta1+theta2),2)) + (pow(sin(theta1-theta2),2)/pow(sin(theta1+theta2),2)) );
						Pt = 1/2.0f*( ((sin(2.0f*theta1)*sin(2.0f*theta2))/(pow(sin(theta1+theta2),2)*pow(cos(theta1-theta2),2))) + ((sin(2.0f*theta1)*sin(2.0f*theta2))/pow(sin(theta1+theta2),2)) );
					}

				}


				// normalize Pr and Pt
				temp_norm = Pr + Pt;
				Pr = Pr/temp_norm;
				Pt = Pt/temp_norm;


				if(ranecu(seed) < Pr)				// reflection
				{
					trans_dir_cos(dcos, normal, theta1, theta2, 0, info);  // compute reflected directional cosines


					// condition to check that reflected vector is within 1.57 radians from original normal
					angle_oldN_R = dot_product(old_normal, dcos);
					angle_oldN_R = acosf(angle_oldN_R);


					if (angle_oldN_R > 1.57f) // > 1.57 radians, reperturb the normal
					{
						reperturb_ctr++;

						if(reperturb_ctr < 4)		// maximum 3 times reperturb
							goto reperturb;
						else				// calculate using smooth surface normal (old_normal)
						{
							normal[0] = old_normal[0];
							normal[1] = old_normal[1];
							normal[2] = old_normal[2];

							dcos[0] = old_dcos[0];
							dcos[1] = old_dcos[1];
							dcos[2] = old_dcos[2];

							dcos_temp[0] = -dcos[0];
							dcos_temp[1] = -dcos[1];
							dcos_temp[2] = -dcos[2];

							reperturb_ctr = 0;

							if(oldN_Rctr < 25)	// max resample 25 times (25*3 reperturb = 75 times)
							{
								oldN_Rctr++;
								goto no_perturbation;
							}
							else			// terminate photon
							{
								num_theta1++;	// increment the counter for 'theta1' - # photons terminated due to incidence angle > 1.57 or < 0 radian
								flag_abs = 1;
								oldN_Rctr = 0;
								goto baexit;
							}
						}
					}

					if (flag_call_transmit == 0)		// reflects between columns, calculate distance using transmit()
					{
						flagCCT = 0;

						trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

						if (trans_flag == 1)		// photon terminated
							flag_abs = 1;
						else if (trans_flag == 0)
							goto prpt;				
					}
				}
				else						// transmission
				{
					trans_dir_cos(dcos, normal, theta1, theta2, 1, info);	// compute transmitted directional cosines

					if (flag_call_transmit == 1)		// photon exits current column
					{
						flag_call_transmit = 0;
						flagCCT = 0;

						trans_flag = transmit(pos, dcos, normal, seed, xdetector, ydetector, H, top_absfrac, beta, d_min, pixelsize, lbound_x, lbound_y, ubound_x, ubound_y, myimage, info, d_max, sensorRefl, ydim, flagCCT, h_num_detected_prim);

						if (trans_flag == 1)		// photon terminated
							flag_abs = 1;
						else if (trans_flag == 0)	// hits a column
							goto prpt;		// check again to see if it gets reflected or transmitted
					}			
				}
			} // else prpt ends
		}	// else ends

	baexit:
	   return flag_abs;

	}	

#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Photon gets transmitted, calculate the new position where it hits next column or boundary
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
__device__ int transmit(float3 *pos, float3 *dcos, float3 *normal, int2* seed, float xdetector, float ydetector, float H, 
float top_absfrac, float beta, float d_min, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y,
unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float d_max, 
float sensorRefl, int flagCCT, float *Xc, float *Yc, struct histStruct * dPhotonHist, int boolToCollectData)
{
	float d_nextCol = 0.0f;		// distance to next column	
	float d_top = 0.0f;		// distance to top surface
	float d_bottom = 0.0f;		// distance to bottom surface (sensor plane)
	int photon_exit = 0;		// flag - photon terminates (lost/detected/absorbed) (yes=1; no=0)
	int reflbtm = 0;		// flag - photon terminated during call to refl_bottom() (yes=1; no=0)
	int newnormalctr=0;		// counter - if angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate normal; max. 100 times)
	int newnormalctr2=0;
	float newangle = 0.0f;
	float cos_newangle = 0.0f;
	float3 temp_pos = {0.0f};
	float3 temp_dcos = {0.0f};
	float rr_rnd = 0.0f, theta_rnd = 0.0f;
	float tmp_deno = 0.0f;
	int iii = 0, jjj = 0;

	temp_pos.x = pos->x;
	temp_pos.y = pos->y;
	temp_pos.z = pos->z;

	temp_dcos.x = -dcos->x;		// inverted directional cosines - used to calculate angle between incident and normal vectors
	temp_dcos.y = -dcos->y;
	temp_dcos.z = -dcos->z;

	/*if(mytid == 6 || mytid == 7)
	{
			printf("tid: %d in transmit\n", mytid);
	}*/
	
	if(flagCCT == 1)	// columnar crosstalk occurs
	{
		// photon moves to adjacent column with no change in its dir. cosines. d_nextcol = 0. adjacent column has random orientation.
		newnormal1:
			normal->x = temp_dcos.x;		// invert incident dir. cosines
			normal->y = temp_dcos.y;
			normal->z = temp_dcos.z;

			RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

			tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal->z = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
			normal->x = normal->x/tmp_deno;		
			normal->y = normal->y/tmp_deno;

			// perturb the normal according to Beta
			if(beta > 0.0f)
				RoughSurface(normal, seed, beta);

			// find the angle between normal and -dcos
			cos_newangle = dot_product(&temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// new normal within +- 1.57 radians from inverted dcos
			{						// keep looping until 'newangle' within 1.57 radians (max. 100 iterations)
				if(newnormalctr < 100)
				{
					newnormalctr++;
					goto newnormal1;			
				}
				else 					// else terminate photon
				{
					atomicAdd(&num_theta1,1);	// increment counter - # photons terminated due to incidence angle > 1.57 or < 0 radian
					photon_exit = 1;
					newnormalctr = 0;
					goto exitnow;
				}						
			}
		
			photon_exit = 0;	// photon still alive
	}
	else		// no columnar crosstalk
	{
		// sample distance uniformly between d_min and d_max to next column
		d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

		// new position of the photon. 
		pos->x = temp_pos.x + dcos->x * d_nextCol;
		pos->y = temp_pos.y + dcos->y * d_nextCol;
		pos->z = temp_pos.z + dcos->z * d_nextCol;
	
		if(mytid < dev_numPhotonHist && boolToCollectData == 1)
		{
			dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
			dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
			dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
			dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
			dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
			dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
			dPhotonHist[mytid].histCounter ++;
		}
		
		/*if(mytid == 6 || mytid == 7)
		{
			printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
		}*/
			
		d_top = ((H/2.0f) - temp_pos.z)/dcos->z;
		d_bottom  = ((-H/2.0f) - temp_pos.z)/dcos->z;

		// new position within detector boundaries? - if false, photon LOST
		if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
		{
			atomicAdd(&num_lost, 1);	// increment # photons lost
			
			if(mytid < dev_numPhotonHist && boolToCollectData == 1)
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
				dPhotonHist[mytid].histCounter ++;
			}
			
			/*if(mytid == 6 || mytid == 7)
			{
				printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
			}*/
			
			photon_exit = 1;
			goto exitnow;
		}

		if ( (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
		{
			if( (d_top < d_nextCol) && (d_top > epsilon) )	// d_top < d_nextCol: photon reflects from top surface
			{
				pos->x = temp_pos.x + dcos->x * d_top;
				pos->y = temp_pos.y + dcos->y * d_top;
				pos->z = H/2.0f;

				if(mytid < dev_numPhotonHist && boolToCollectData == 1)
				{
					dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
					dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
					dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
					dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
					dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
					dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
					dPhotonHist[mytid].histCounter ++;
				}
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
			
				atomicAdd(&photon_distance, d_top);	// add distance to mean free path of photon
				photon_exit = 0;			
			}
			else if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )	// d_bottom < d_nextCol, photon hits sensor plane and is reflected/detected
			{
				pos->x = temp_pos.x + dcos->x * d_bottom;
				pos->y = temp_pos.y + dcos->y * d_bottom;
				pos->z = -H/2.0f;

				if(mytid < dev_numPhotonHist && boolToCollectData == 1)
				{
					dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
					dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
					dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
					dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
					dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
					dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
					dPhotonHist[mytid].histCounter ++;
				}
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
			
				atomicAdd(&photon_distance, d_bottom);

				// non-ideal sensor - reflects back sensorRefl% of photons into the current column; detects rest
				if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
				{
					photon_exit = 0;

					// normal pointing +z-direction
					normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

					// obtain reflected dcos from the bottom (specular reflection; 
					// bottom surface is smooth, so normal is not perturbed)
					trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, mytid, info);	// reflection only, so 'refl_theta,trans_theta' = 0

					// sample new distance and place new column
					reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max, mytid, Xc, Yc, dPhotonHist, boolToCollectData);

					if(reflbtm == 1)	// photon terminated in refl_bottom()
					{
						photon_exit = 1;
						goto exitnow;
					}

					// hits top surface after reflecting back
					if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
					{
						goto mytopsurface;
					}
				}
				else		// does not reflect back into column; detected
				{
					photon_exit = 1;	
					atomicAdd(&num_detect, 1);

					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 1;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
			
					iii = floor((pos->x-lbound_x)/pixelsize);	// x,y pixel number in the sensor plane of the detector
					jjj = floor((pos->y-lbound_y)/pixelsize);

					// if photon is detected within lower and upper bounds of point response function: accumulate the signal contribution
					if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
					{	
						unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
						atomicAdd(&current_img[jjj],1);
					}

					atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);	// start array from 0.str_histnum starts from 1
				
					goto exitnow;
				}	
			}
			else	// terminate photon
			{
				atomicAdd(&num_lost, 1);

				if(mytid < dev_numPhotonHist && boolToCollectData == 1)
				{
					dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
					dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
					dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
					dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
					dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
					dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
					dPhotonHist[mytid].histCounter ++;
				}
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
			
				photon_exit = 1;
				goto exitnow;
			}
		}
		else
			atomicAdd(&photon_distance, d_nextCol);		// add distance to mean free path


			// sample new normal to determine orientation of new column.
		newnormal:
			normal->x = temp_dcos.x;		// invert dcos of incident vector
			normal->y = temp_dcos.y;
			normal->z = temp_dcos.z;

			RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

			tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal->z = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
			normal->x = normal->x/tmp_deno;		
			normal->y = normal->y/tmp_deno;

			// perturb the normal according to Beta
			if(beta > 0.0f)
				RoughSurface(normal, seed, beta);

			// angle between normal and -dcos
			cos_newangle = dot_product(&temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 1.57 radians from inverted dcos
			{						// keep looping until 'newangle' within 1.57 radians (max. 100 times)
				if(newnormalctr < 100)
				{
					newnormalctr++;
					goto newnormal;			
				}
				else // kill it
				{
					atomicAdd(&num_theta1,1);
					photon_exit = 1;
					newnormalctr = 0;
					goto exitnow;
				}
						
			}
	
			// check if the photon enters another column or get lost (hit detector side)/ reflected (at top surface)/ detected (at bottom surface)
	
			// hits side of detector?
			if ( (fabs(pos->x-0.0f) < epsilon) || (fabs(pos->x-xdetector) < epsilon) || (fabs(pos->y-0.0f) < epsilon) || (fabs(pos->y-ydetector) < epsilon) )		
			{
				atomicAdd(&num_lost, 1);

				if(mytid < dev_numPhotonHist && boolToCollectData == 1)
				{
					dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
					dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
					dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
					dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
					dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
					dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
					dPhotonHist[mytid].histCounter ++;
				}
				
				/*if(mytid == 6 || mytid == 7)
				{
					printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
				}*/
			
				photon_exit = 1;
				goto exitnow;

			}
	

			// hits top surface?
		     mytopsurface:

			if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	// gets reflected or absorbed
			{
				normal->x = 0.0f;
				normal->y = 0.0f;
				normal->z = -1.0f;

				// top surface absorption - using absorption coefficient 'top_absfrac'
				if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )	// photon absorbed
				{
					atomicAdd(&num_abs_top, 1);	// increment # photons absorbed at top counter
	
					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 3;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
			
					photon_exit = 1;
					goto exitnow;
				}
				else	// not absorbed at top
				{
					// assign new directional cosines (isotropic surface)
					dcos->z = -fabs((ranecu(seed) * 2.0f) - 1.0f);
					rr_rnd = sqrt(1.0f - dcos->z*dcos->z);
					theta_rnd = ranecu(seed)*twopipen;	
	
					dcos->x=rr_rnd*cos(theta_rnd);
					dcos->y=rr_rnd*sin(theta_rnd);

					temp_pos.x = pos->x;
					temp_pos.y = pos->y;
					temp_pos.z = pos->z;

					temp_dcos.x = -dcos->x;
					temp_dcos.y = -dcos->y;
					temp_dcos.z = -dcos->z;

					// sample distance uniformly between d_min and d_max to next column
					d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

					// distance to bottom surface: if d_bottom < d_nextCol, photon is detected.
					d_bottom  = ((-H/2.0f) - temp_pos.z)/dcos->z;

					// new position of the photon. 
					pos->x = temp_pos.x + dcos->x * d_nextCol;
					pos->y = temp_pos.y + dcos->y * d_nextCol;
					pos->z = temp_pos.z + dcos->z * d_nextCol;
					
					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
			
					// pos is within detector boundaries? - if false, photon lost
					if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
					{
						atomicAdd(&num_lost, 1);	// increment # photons lost
	
						if(mytid < dev_numPhotonHist && boolToCollectData == 1)
						{
							dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
							dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
							dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
							dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
							dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
							dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
							dPhotonHist[mytid].histCounter ++;
						}
						
						/*if(mytid == 6 || mytid == 7)
						{
							printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
						}*/
			
						photon_exit = 1;
						goto exitnow;
					}

					if ( (pos->z < -H/2.0f) || (pos->z > H/2.0f)  )
					{
						if( (d_bottom < d_nextCol) && (d_bottom > epsilon) ) // dist. to bottom < dist. to next column
						{
							pos->x = temp_pos.x + dcos->x * d_bottom;
							pos->y = temp_pos.y + dcos->y * d_bottom;
							pos->z = -H/2.0f;
			
							if(mytid < dev_numPhotonHist && boolToCollectData == 1)
							{
								dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
								dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
								dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
								dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
								dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
								dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 0;
								dPhotonHist[mytid].histCounter ++;
							}
							
							/*if(mytid == 6 || mytid == 7)
							{
								printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
							}*/
			
							atomicAdd(&photon_distance, d_bottom);

							// non-ideal sensor - reflects back sensorRefl% of photons into the current column; detects rest
							if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
							{
								photon_exit = 0;

								// normal pointing (0,0,1)
								normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

								// obtain reflected dcos from the bottom (specular reflection; 
								// bottom surface is smooth, do not perturb the normal)
								trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, mytid, info);   // reflection only so 'refl_theta,trans_theta' = 0

								// sample new distance and place new column there
								reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max, mytid, Xc, Yc, dPhotonHist, boolToCollectData);

								if(reflbtm == 1)
								{
									photon_exit = 1;
									goto exitnow;
								}

								// hits top surface after reflecting back
								if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
								{
									goto mytopsurface;
								}
							}
							else	// not reflected back into column; detected
							{
								photon_exit = 1;
								atomicAdd(&num_detect, 1);

								if(mytid < dev_numPhotonHist && boolToCollectData == 1)
								{
									dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
									dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
									dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
									dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
									dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
									dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 1;
									dPhotonHist[mytid].histCounter ++;
								}
								
								/*if(mytid == 6 || mytid == 7)
								{
									printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
								}*/
				
								iii = floor((pos->x-lbound_x)/pixelsize);	// x,y pixel number
								jjj = floor((pos->y-lbound_y)/pixelsize);

								// if photon gets detected within lower and upper bounds: accumulate signal contribution
								if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
								{	
									unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
									atomicAdd(&current_img[jjj],1);
								}

								atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);// increment # detected per primary

								goto exitnow;	
							}	
						}
						else	// terminate photon
						{
							atomicAdd(&num_lost, 1);

							if(mytid < dev_numPhotonHist && boolToCollectData == 1)
							{
								dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
								dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
								dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
								dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
								dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
								dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
								dPhotonHist[mytid].histCounter ++;
							}
							
							/*if(mytid == 6 || mytid == 7)
							{
								printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
							}*/
			
							photon_exit = 1;
							goto exitnow;
						}
					}
					else	// photon lies between top and bottom surfaces of detector
						atomicAdd(&photon_distance, d_nextCol);		// add distance travelled to global variable

					// sample new normal to determine orientation of new column.
			  newnormal_TOP:
					normal->x = temp_dcos.x;		// invert dcos of incident vector
					normal->y = temp_dcos.y;
					normal->z = temp_dcos.z;

					RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

					tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

					// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
					normal->z = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
					normal->x = normal->x/tmp_deno;		
					normal->y = normal->y/tmp_deno;

					// perturb the normal according to beta
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);

					// angle between normal and -dcos
					cos_newangle = dot_product(&temp_dcos, normal);
					newangle = acosf(cos_newangle);

					if ( (newangle < 0.0f) || (newangle > 1.57f) )	// new normal within +- 1.57 radians from inverted dcos
					{						// keep looping until 'newangle' within 1.57 radians (max. 100 times)
						if(newnormalctr2 < 100)			
						{
							newnormalctr2++;
							goto newnormal_TOP;					
						}
						else 					// terminate photon
						{
							atomicAdd(&num_theta1,1);
							photon_exit = 1;
							newnormalctr2 = 0;
							goto exitnow;
						}
					}
					photon_exit = 0;
				}
			}	// hit top ends
	

			// hit bottom? 
			if ( fabs(pos->z - (-H/2.0f)) < epsilon )	// gets detected
			{
				// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
				if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
				{
					photon_exit = 0;

					// normal pointing (0,0,1)
					normal->x = 0.0f; normal->y = 0.0f; normal->z = 1.0f;

					// obtain reflected dcos from the bottom (specular reflection; 
					// bottom surface is smooth, do not perturb the normal)
					trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, mytid, info);	// reflection only so 'refl_theta,trans_theta' = 0

					// sample new distance and place new column
					reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max, mytid, Xc, Yc, dPhotonHist, boolToCollectData);

					if(reflbtm == 1)
					{
						photon_exit = 1;
						goto exitnow;
					}

					// hits top surface after reflecting back?
					if ( (fabs(pos->z - (H/2.0f)) < epsilon) && (dcos->z > 0.0f) )	
					{
						goto mytopsurface;
					}
				}
				else	// detected
				{
					atomicAdd(&num_detect, 1);

					if(mytid < dev_numPhotonHist && boolToCollectData == 1)
					{
						dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
						dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
						dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
						dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
						dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
						dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 1;
						dPhotonHist[mytid].histCounter ++;
					}
					
					/*if(mytid == 6 || mytid == 7)
					{
						printf("transmit - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
					}*/
			
					photon_exit = 1;

					iii = floor((pos->x-lbound_x)/pixelsize);	// x,y pixel number in the sensor plane of the detector
					jjj = floor((pos->y-lbound_y)/pixelsize);

					// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
					if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
					 {	
						unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + iii * pitch);
						atomicAdd(&current_img[jjj],1);
					 }

					atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);	// increment # detected per primary
	
					goto exitnow;
				}
			}
		} // else flagCCT ends
	
	exitnow:
	 return photon_exit;	
	}
	
#else	// C code

	int transmit(float *pos, float *dcos, float *normal, int* seed, float xdetector, float ydetector, float H, float top_absfrac, float beta, float d_min, int pixelsize, 
	float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, float d_max, float sensorRefl, int ydim, int flagCCT, 		int *h_num_detected_prim)
	{
		float d_nextCol = 0.0f;		// distance to next column	
		float d_top = 0.0f;		// distance to top surface
		float d_bottom = 0.0f;		// distance to bottom surface
		int photon_exit = 0;		// flag - photon terminates (lost/detected/absorbed) (yes=1; no=0)
		int reflbtm = 0;		// flag - photon terminated during call to refl_bottom() (yes=1; no=0)
		int newnormalctr = 0;		// counter - if angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate normal; max. 100 times)
		int newnormalctr2 = 0;
		float newangle = 0.0f;
		float cos_newangle = 0.0f;
		float temp_pos[3] = {0.0f};
		float temp_dcos[3] = {0.0f};
		float rr_rnd = 0.0f, theta_rnd = 0.0f;
		float tmp_deno = 0.0f;
		int iii = 0, jjj = 0;


		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		temp_dcos[0] = -dcos[0];	// inverted directional cosine - used to calculate angle between incident vector and normal
		temp_dcos[1] = -dcos[1];
		temp_dcos[2] = -dcos[2];

		if(flagCCT == 1)	// columnar crosstalk occurs
		{
			// cross over to adjacent column. no change in dir. cosines. d_nextcol = 0. adjacent column has random orientation.
			newnormal1:
				normal[0] = temp_dcos[0];		// invert dcos of incident vector
				normal[1] = temp_dcos[1];
				normal[2] = temp_dcos[2];

				RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

				tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal[2] = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
				normal[0] = normal[0]/tmp_deno;		
				normal[1] = normal[1]/tmp_deno;

				// perturb the normal according to Beta
				if(beta > 0.0f)
					RoughSurface(normal, seed, beta);

				// find the angle between Normal and -Dcos
				cos_newangle = dot_product(temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 1.57 radians from inverted dcos
				 {						// keep looping until 'newangle' within 1.57 radians (max. 100 times)
					if(newnormalctr < 100)
					{
						newnormalctr++;
						goto newnormal1;			
					}
					else 					// terminate photon
					{
						num_theta1++;
						photon_exit = 1;
						newnormalctr = 0;
						goto exitnow;
					}
						
				 }

				photon_exit = 0;	// photon still alive
		}
		else		// no columnar crosstalk occurs
		{

			// sample distance uniformly between d_min and d_max to next column
			d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

			// new position of the photon. 
			pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
			pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
			pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

			d_top = ((H/2.0f) - temp_pos[2])/dcos[2];
			d_bottom  = ((-H/2.0f) - temp_pos[2])/dcos[2];

			// condition to check that pos is within detector boundaries - if true, photon LOST
			if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
			{
				num_lost++;
				photon_exit = 1;
				goto exitnow;
			}

			if ( (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
				{
					if( (d_top < d_nextCol) && (d_top > epsilon) )	// d_top < d_nextCol: photon reflects from the top surface
					{
						pos[0] = temp_pos[0] + dcos[0] * d_top;
						pos[1] = temp_pos[1] + dcos[1] * d_top;
						pos[2] = H/2.0f;
				
						photon_distance = photon_distance + d_top;
						photon_exit = 0;			
					}
					else if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )	// d_bottom < d_nextCol: photon detected.
					{
						pos[0] = temp_pos[0] + dcos[0] * d_bottom;
						pos[1] = temp_pos[1] + dcos[1] * d_bottom;
						pos[2] = -H/2.0f;

						photon_distance = photon_distance + d_bottom;

						// non-ideal sensor - reflects back sensorRefl% of photons into the current column; detects rest
						if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
						{
		
							photon_exit = 0;

							// normal pointing +z-direction
							normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

							// obtain reflected dcos from the bottom (specular reflection; 
							// bottom surface is smooth, do not perturb the normal)
							trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so 'refl_theta,trans_theta' = 0

							// sample new distance and place new column there
							reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

							if(reflbtm == 1)
							{
								photon_exit = 1;	// photon terminated in refl_bottom()
								goto exitnow;
							}


							// hits top surface after reflecting back?
							if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
							{
								goto mytopsurface;
							}

						}
						else	// not reflected back into column; photon detected
						{
							photon_exit = 1;	
							num_detect++;

							iii = floor((pos[0]-lbound_x)/pixelsize);	// x,y pixel number in the sensor plane of the detector
							jjj = floor((pos[1]-lbound_y)/pixelsize);

							// if the photon gets detected within lower and upper bounds: accumulate the signal contribution
							if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
							 {	
								outputimage_.newimageopt[iii][jjj]++;
							 }

							h_num_detected_prim[info.str_histnum]++;	// increment # detected per primary
				
							goto exitnow;
						}	
					}
					else	// terminate photon
					{
						num_lost++;
						photon_exit = 1;
						goto exitnow;
					}
				}
			else		// photon lies between top and bottom surfaces of detector
				photon_distance = photon_distance + d_nextCol;		// add distance to mean free path of photon


			// sample new normal to determine orientation of new column.
		newnormal:
			normal[0] = temp_dcos[0];		// invert dcos of incident vector
			normal[1] = temp_dcos[1];
			normal[2] = temp_dcos[2];

			RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

			tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal[2] = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
			normal[0] = normal[0]/tmp_deno;		
			normal[1] = normal[1]/tmp_deno;

			// perturb the normal according to Beta
			if(beta > 0.0f)
				RoughSurface(normal, seed, beta);

			// angle between normal and -dcos
			cos_newangle = dot_product(temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 1.57 radians from inverted dcos
			 {						// keep looping until 'newangle' within 1.57 radians (max. 100 iterations)
				if(newnormalctr < 100)
				{
					newnormalctr++;
					goto newnormal;			
				}
				else 					// terminate photon
				{
					num_theta1++;
					photon_exit = 1;
					newnormalctr = 0;
					goto exitnow;
				}
						
			 }
	
			// check if the photon enters another column or gets lost (hit detector side)/ reflected (detector top)/ detected (detector bottom)
	
			// hits side of detector?
			if ( (fabs(pos[0]-0.0f) < epsilon) || (fabs(pos[0]-xdetector) < epsilon) || (fabs(pos[1]-0.0f) < epsilon) || (fabs(pos[1]-ydetector) < epsilon) )		
			{
				num_lost++;
				photon_exit = 1;
				goto exitnow;

			}


			// hit top?
		     mytopsurface:

			if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	// gets reflected or absorbed
			{
				normal[0] = 0.0f;
				normal[1] = 0.0f;
				normal[2] = -1.0f;

				if ( (top_absfrac > 0.0f) && (ranecu(seed) < top_absfrac) )	// photon gets absorbed at top surface
				{
					num_abs_top++;
					photon_exit = 1;
					goto exitnow;
				}
				else	// photon reflected from top
				{
					// assign new directional cosines (isotropic top surface)
					dcos[2] = -fabs((ranecu(seed) * 2.0f) - 1.0f);
					rr_rnd = sqrt(1.0f - dcos[2]*dcos[2]);
					theta_rnd = ranecu(seed)*twopipen;	
	
					dcos[0]=rr_rnd*cos(theta_rnd);
					dcos[1]=rr_rnd*sin(theta_rnd);

					temp_pos[0] = pos[0];
					temp_pos[1] = pos[1];
					temp_pos[2] = pos[2];

					temp_dcos[0] = -dcos[0];
					temp_dcos[1] = -dcos[1];
					temp_dcos[2] = -dcos[2];

					// sample distance uniformly between d_min and d_max to next column
					d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

					// distance to bottom surface: if d_bottom < d_nextCol, photon detected.
					d_bottom  = ((-H/2.0f) - temp_pos[2])/dcos[2];

					// compute the new position of the photon. 
					pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
					pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
					pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

					// check new position is within detector boundaries? - if false, photon lost
					if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
					{
						num_lost++;
						photon_exit = 1;
						goto exitnow;
					}

					if ( (pos[2] < -H/2.0f) || (pos[2] > H/2.0f)  )
						{
							if( (d_bottom < d_nextCol) && (d_bottom > epsilon) )
							{
								pos[0] = temp_pos[0] + dcos[0] * d_bottom;
								pos[1] = temp_pos[1] + dcos[1] * d_bottom;
								pos[2] = -H/2.0f;

								photon_distance = photon_distance + d_bottom;

								// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
								if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
								{

									photon_exit = 0;		
			
									// normal pointing (0,0,1)
									normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

									// obtain reflected dcos from the bottom (specular reflection; 
									// bottom surface is smooth, do not perturb the normal)
									trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so 'refl_theta,trans_theta' = 0

									// sample new distance and place new column
									reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

									if(reflbtm == 1)
									{
										photon_exit = 1;
										goto exitnow;
									}

									// hits top surface after reflecting back?
									if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
									{
										goto mytopsurface;
									}

								}
								else		// not reflected back into column; detected
								{
										photon_exit = 1;
										num_detect++;

										iii = floor((pos[0]-lbound_x)/pixelsize);	// x,y pixel number
										jjj = floor((pos[1]-lbound_y)/pixelsize);

										// photon gets detected within lower and upper bounds: accumulate signal contribution
										if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
										 {	
											outputimage_.newimageopt[iii][jjj]++;
										 }

										h_num_detected_prim[info.str_histnum]++;	// increment # detected per primary

										goto exitnow;	
								}	
							}
							else	// terminate photon
							{
								num_lost++;
								photon_exit = 1;
								goto exitnow;
							}
						}
					else
						photon_distance = photon_distance + d_nextCol;		// add distance to mean free path of photon

					// sample new normal to determine orientation of new column.
			  newnormal_TOP:
					normal[0] = temp_dcos[0];		// invert dcos of incident vector
					normal[1] = temp_dcos[1];
					normal[2] = temp_dcos[2];

					RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

					tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

					// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
					normal[2] = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
					normal[0] = normal[0]/tmp_deno;		
					normal[1] = normal[1]/tmp_deno;

					// perturb the normal according to Beta
					if(beta > 0.0f)
						RoughSurface(normal, seed, beta);

					// angle between normal and -dcos
					cos_newangle = dot_product(temp_dcos, normal);
					newangle = acosf(cos_newangle);

					if ( (newangle < 0.0f) || (newangle > 1.57f) )	// new normal within +- 1.57 radians from inverted dcos
					 {						// keep looping until 'newangle' within 1.57 radians (max. 100 iterations)
						if(newnormalctr2 < 100)	
						{
							newnormalctr2++;
							goto newnormal_TOP;			
						}
						else 					// terminate photon
						{
							num_theta1++;
							photon_exit = 1;
							newnormalctr2 = 0;
							goto exitnow;
						}
				
					 }
		
					photon_exit = 0;	// photon still alive
				}
			}	// hit top ends
	

			// hit bottom? 
			if ( fabs(pos[2] - (-H/2.0f)) < epsilon )	// gets detected
			{
				// non-ideal sensor - reflects back sensorRefl% of photons into the current column; absorbs rest
				if(ranecu(seed) < sensorRefl)	// reflect back - specular (mirror) reflection
				{
					photon_exit = 0;		
	
					// normal pointing (0,0,1)
					normal[0] = 0.0f; normal[1] = 0.0f; normal[2] = 1.0f;

					// obtain reflected dcos from the bottom (specular reflection; 
					// bottom surface is smooth, do not perturb the normal)
					trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, info);	// reflection only so 'refl_theta,trans_theta' = 0

					// sample new distance and place new column
					reflbtm = refl_bottom(pos, dcos, normal, xdetector, ydetector, seed, beta, d_min, H, d_max);

					if(reflbtm == 1)
					{
						photon_exit = 1;
						goto exitnow;
					}

					// hits top surface after reflecting back?
					if ( (fabs(pos[2] - (H/2.0f)) < epsilon) && (dcos[2] > 0.0f) )	
					{
						goto mytopsurface;
					}

				}
				else		// detected
				{
					num_detect++;
					photon_exit = 1;

					iii = floor((pos[0]-lbound_x)/pixelsize);	// x,y pixel number in sensor plane of detector
					jjj = floor((pos[1]-lbound_y)/pixelsize);

					// photon gets detected within lower and upper bounds: accumulate the signal contribution
					if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
					 {	
						outputimage_.newimageopt[iii][jjj]++;
					 }

					h_num_detected_prim[info.str_histnum]++;	// increment # detected per primary
	
					goto exitnow;
				}
			}

		} // else flagCCT ends

	exitnow:
	 return photon_exit;	
	}
	
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// photon reflects from sensor_plane or bottom surface, when in between columns. 
// Obtains the next column where it hits.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
__device__ int refl_bottom(float3 *pos, float3 *dcos, float3 *normal, float xdetector, float ydetector, int2* seed, 
float beta, float d_min, float H, float d_max, int mytid, float *Xc, float *Yc, struct histStruct * dPhotonHist, int boolToCollectData)
{
	float d_nextCol=0.0f;		// distance to next column
	float d_top=0.0f;		// distance to top surface of detector
	int pexit=0;			// flag - photon terminated in this function call (exited: yes=1; no=0)
	int newnormalctr=0;		// counter - angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate max. 100 times)
	float tmp_deno=0.0f, cos_newangle=0.0f, newangle=0.0f;
	float3 temp_pos, temp_dcos;

	temp_pos.x = pos->x;
	temp_pos.y = pos->y;
	temp_pos.z = pos->z;

	temp_dcos.x = -dcos->x;		// inverted directional cosines - used to calculate angle between incident and normal vectors
	temp_dcos.y = -dcos->y;
	temp_dcos.z = -dcos->z;

	// sample distance uniformly between d_min and d_max to next column
	d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

	// distance to top surface
	d_top  = ((H/2.0f) - temp_pos.z)/dcos->z;

	// new position of the photon - specular reflection 
	pos->x = temp_pos.x + dcos->x * d_nextCol;
	pos->y = temp_pos.y + dcos->y * d_nextCol;
	pos->z = temp_pos.z + dcos->z * d_nextCol;

	/*if(mytid == 6 || mytid == 7)
	{
			printf("tid: %d in refl_bottom\n", mytid);
	}*/
	
	// is new position is within detector boundaries? if false, photon lost
	if ( (pos->x < epsilon) || (pos->x > xdetector) || (pos->y < epsilon) || (pos->y > ydetector) )
	{
		atomicAdd(&num_lost, 1);	// increment # photons lost

		if(mytid < dev_numPhotonHist && boolToCollectData == 1)
		{
			dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
			dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
			dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
			dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
			dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
			dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
			dPhotonHist[mytid].histCounter ++;
		}
		/*if(mytid == 6 || mytid == 7)
		{
			printf("refl_bottom - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
		}*/
			
		pexit = 1;			// set the flag
		goto myexit;
	}

	if ( (pos->z > H/2.0f)  )		// if photon's new z-position is above the top surface of detector
	{
		if( (d_top < d_nextCol) && (d_top > epsilon) )		// if distance to top < dist. to next column - photon will hit top surface
		{
			pos->x = temp_pos.x + dcos->x * d_top;
			pos->y = temp_pos.y + dcos->y * d_top;
			pos->z = H/2.0f;
				
			atomicAdd(&photon_distance, d_top);		// add this distance to mean free path of photon
			pexit = 0;			
		}
		else			// terminate photon
		{
			atomicAdd(&num_lost, 1);

			if(mytid < dev_numPhotonHist && boolToCollectData == 1)
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 2;
				dPhotonHist[mytid].histCounter ++;
			}
			
			/*if(mytid == 6 || mytid == 7)
			{
				printf("refl_bottom - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
			}*/
			
			pexit = 1;
			goto myexit;
		}
	}
	else					// photon z coordinate is less than H/2
	{
		atomicAdd(&photon_distance, d_nextCol);		// add distance to mean free path of photon

		// sample new normal to determine orientation of new column.
		 mynewnormal:
			normal->x = temp_dcos.x;		// invert incident directional cosines
			normal->y = temp_dcos.y;
			normal->z = temp_dcos.z;

			RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

			tmp_deno = sqrt(normal->x*normal->x + normal->y*normal->y);

			// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
			normal->z = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
			normal->x = normal->x/tmp_deno;		
			normal->y = normal->y/tmp_deno;

			// perturb the normal according to Beta
			RoughSurface(normal, seed, beta);

			// angle between normal and -dcos
			cos_newangle = dot_product(&temp_dcos, normal);
			newangle = acosf(cos_newangle);

			if ( (newangle < 0.0f) || (newangle > 1.57f) )	// new rough normal is within +- 1.57 radians from inverted dcos
			{
				if(newnormalctr < 100)	// max 100 times
				{
					newnormalctr++;
					goto mynewnormal;		// keep looping until 'newangle' within 1.57 radians	
				}
				else // kill it
				{
					atomicAdd(&num_theta1,1);   // increment the counter for 'theta1' - # photons terminated due to incidence angle > 1.57 or < 0 radian
					pexit = 1;
					newnormalctr = 0;
					goto myexit;
				}
			}
			pexit = 0;	// photon still alive
	}
	myexit:
		return pexit;
}	
#else	// C code

	int refl_bottom(float *pos, float *dcos, float *normal, float xdetector, float ydetector, int* seed, float beta, float d_min, float H, float d_max)
	{
		float d_nextCol=0.0f;	// distane to next column
		float d_top=0.0f;	// distance to top surface of detector
		int pexit=0;		// flag - photon ternimated in this function call (exited: yes=1; no=0)
		int newnormalctr = 0;	// counter - angle between inverted dir. cosine and rough normal > 1.57 or < 0 radians (recalculate max. 100 times)
		float temp_pos[3], temp_dcos[3];
		float tmp_deno=0.0f, cos_newangle=0.0f, newangle=0.0f;


		temp_pos[0] = pos[0];
		temp_pos[1] = pos[1];
		temp_pos[2] = pos[2];

		temp_dcos[0] = -dcos[0];		// inverted directional cosines - used to calculate angle between incident and normal vectors
		temp_dcos[1] = -dcos[1];
		temp_dcos[2] = -dcos[2];

		// sample distance uniformly between d_min and d_max to next column
		d_nextCol = ranecu(seed) * (d_max - d_min) + d_min;

		// distance to top surface
		d_top  = ((H/2.0f) - temp_pos[2])/dcos[2];

		// compute the new position of the photon. 
		pos[0] = temp_pos[0] + dcos[0] * d_nextCol;
		pos[1] = temp_pos[1] + dcos[1] * d_nextCol;
		pos[2] = temp_pos[2] + dcos[2] * d_nextCol;

		// is new position is within detector boundaries? if false, photon lost
		if ( (pos[0] < epsilon) || (pos[0] > xdetector) || (pos[1] < epsilon) || (pos[1] > ydetector) )
		{
			num_lost++;	// increment # photons lost
			pexit = 1;
			goto myexit;
		}

		if ( (pos[2] > H/2.0f)  )		// photon's new z position is above top surface
		{
				if( (d_top < d_nextCol) && (d_top > epsilon) )		// if distance to top < dist. to next column; photon will hit top surface
				{
					pos[0] = temp_pos[0] + dcos[0] * d_top;
					pos[1] = temp_pos[1] + dcos[1] * d_top;
					pos[2] = H/2.0f;
				
					photon_distance = photon_distance + d_top;	// add this distance to mean free path of photon
					pexit = 0;			
				}
				else			// terminate photon
				{
					num_lost++;
					pexit = 1;
					goto myexit;
				}
		}
		else
		{
			photon_distance = photon_distance + d_nextCol;		// add distance to mean free path of photon

			// sample new normal to determine orientation of new column.
		  	mynewnormal:

				normal[0] = temp_dcos[0];		// invert incident directional cosines
				normal[1] = temp_dcos[1];
				normal[2] = temp_dcos[2];

				RoughSurface(normal, seed, 1.0f);	// beta = 1.0 to get new normal within +-1.57 radians of inverted dcos.

				tmp_deno = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

				// make z component of normal equal to 0.0f and renormalize (because we want to rotate the normal only in the x-y plane)
				normal[2] = 0.0f;			// normal_z of a cylinder is zero (no tilt assumed)
				normal[0] = normal[0]/tmp_deno;		
				normal[1] = normal[1]/tmp_deno;

				// perturb the normal according to Beta
				RoughSurface(normal, seed, beta);

				// angle between Normal and -Dcos
				cos_newangle = dot_product(temp_dcos, normal);
				newangle = acosf(cos_newangle);

				if ( (newangle < 0.0f) || (newangle > 1.57f) )	// check if new normal is within +- 1.57 radians from inverted dcos
				 {
					if(newnormalctr < 100)	// max 100 times
					{
						newnormalctr++;
						goto mynewnormal;			// keep looping until 'newangle' within 1.57 radians	
					}
					else 			// terminate photon
					{
						num_theta1++;	// increment the counter for 'theta1' - # photons terminated due to incidence angle > 1.57 or < 0 radian
						pexit = 1;
						newnormalctr = 0;
						goto myexit;
					}		
				 }
				pexit = 0;	// photon still alive
		}

	myexit:
	 return pexit;
	}
	
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// calculate dot product of two vectors
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

	__device__ inline float dot_product(float3 *aa, float3 *b)
	{
		float result = 0.0f;

		result = aa->x*b->x + aa->y*b->y + aa->z*b->z;

	  return result;
	}

#else	// C code

	float dot_product(float *aa, float *b)
	{
		float result = 0.0f;

		result = aa[0]*b[0] + aa[1]*b[1] + aa[2]*b[2];

	  return result;
	}
	
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// compute directional cosines of transmitted vector
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

	__device__ inline void trans_dir_cos(float3 *dcos, float3 *normal, float refl_theta, float trans_theta, int flag_ref, int mytid, struct start_info *info)
	{
		float cos_angle = 0.0f;
		float norm = 0.0f;
		float3 dcos_temp = {0.0f};

		dcos_temp.x = -dcos->x;		// inverted directional cosines - used to calculate angle between incident and normal vectors
		dcos_temp.y = -dcos->y;
		dcos_temp.z = -dcos->z;
	
		cos_angle = dot_product(&dcos_temp,normal);	// cosine of angle between inverted incident dir. cosine and normal

		if (flag_ref == 0)				// reflection
		{
			dcos->x = 2.0f*cos_angle*normal->x + dcos->x;  // specular reflection
			dcos->y = 2.0f*cos_angle*normal->y + dcos->y;
			dcos->z = 2.0f*cos_angle*normal->z + dcos->z;
		}
		else if (flag_ref == 1)				// transmission	
		{
			 dcos->x= -normal->x*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->x+(cos_angle*normal->x));
			 dcos->y= -normal->y*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->y+(cos_angle*normal->y));
			 dcos->z= -normal->z*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos->z+(cos_angle*normal->z));
		}

		// normalize
		norm = sqrt(dcos->x*dcos->x + dcos->y*dcos->y + dcos->z*dcos->z);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
		 {
			dcos->x = dcos->x/norm;
			dcos->y = dcos->y/norm;
			dcos->z = dcos->z/norm;
		 } 

	return;	
	}

#else	// C code

	void trans_dir_cos(float *dcos, float *normal, float refl_theta, float trans_theta, int flag_ref, struct start_info info)
	{
		float cos_angle = 0.0f;
		float norm = 0.0f;
		float dcos_temp[3] = {0.0f};

		dcos_temp[0] = -dcos[0];	// inverted directional cosines - used to calculate angle between incident and normal vectors
		dcos_temp[1] = -dcos[1];
		dcos_temp[2] = -dcos[2];
	
		cos_angle = dot_product(dcos_temp,normal);	// cosine of angle between inverted incident dir. cosine and normal

		if (flag_ref == 0)				// reflection
		{
				dcos[0] = 2.0f*cos_angle*normal[0] + dcos[0];  // specular reflection
				dcos[1] = 2.0f*cos_angle*normal[1] + dcos[1];
				dcos[2] = 2.0f*cos_angle*normal[2] + dcos[2];
		}
		else if (flag_ref == 1)				// transmission	
		{
			 dcos[0]= -normal[0]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[0]+(cos_angle*normal[0]));
			 dcos[1]= -normal[1]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[1]+(cos_angle*normal[1]));
			 dcos[2]= -normal[2]*cos(trans_theta)-(sin(trans_theta)/sin(refl_theta))*(dcos[2]+(cos_angle*normal[2]));
		}

		// normalize
		norm = sqrt(dcos[0]*dcos[0] + dcos[1]*dcos[1] + dcos[2]*dcos[2]);

		if ((norm < (1.0f - epsilon)) || (norm > (1.0f + epsilon)))
		 {
			dcos[0] = dcos[0]/norm;
			dcos[1] = dcos[1]/norm;
			dcos[2] = dcos[2]/norm;
		 } 

	return;	
	}

#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// add roughness to the surface of the column according to roughness coefficient 'beta'
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA

	__device__ inline void RoughSurface(float3 *normal, int2* seed, float beta)
	{

		float theta = 0.0f;
		float status = 0.0f;
		float rr = 0.0f;
		float3 normalpert = {0.0f};
		float3 rough_normal = {0.0f};
		float normalize_base = 0.0f;

		// generate the perturbation vector
		status = ranecu(seed);
		normalpert.z = 2.0f*status - 1.0f;	// random number between (-1,1)
		rr = sqrt(1.0f - status*status);
		status = ranecu(seed);
		theta = status * 2.0f * pi;		// random number between (0,2pi)

		normalpert.x = rr * cos(theta);
		normalpert.y = rr * sin(theta);

		// normalize the perturbed vector
		normalize_base = sqrt( pow(normalpert.x,2) + pow(normalpert.y,2) + pow(normalpert.z,2) );
	
		normalpert.x = normalpert.x/normalize_base;
		normalpert.y = normalpert.y/normalize_base;
		normalpert.z = normalpert.z/normalize_base;

		// rough normal = beta*perturbed + original normal
		rough_normal.x = beta * normalpert.x + normal->x;	
		rough_normal.y = beta * normalpert.y + normal->y;
		rough_normal.z = beta * normalpert.z + normal->z;

		// normalize rough normal
		normalize_base = sqrt( pow(rough_normal.x,2) + pow(rough_normal.y,2) + pow(rough_normal.z,2) );

		normal->x = rough_normal.x/normalize_base; 
		normal->y = rough_normal.y/normalize_base;
		normal->z = rough_normal.z/normalize_base;

	return;
	}
	
#else	// C code

	void RoughSurface(float *normal, int* seed, float beta)
	{
		float theta = 0.0f;
		float status = 0.0f;
		float rr = 0.0f;
		float normalpert[3] = {0.0f};
		float rough_normal[3] = {0.0f};
		float normalize_base = 0.0f;

		// generate the perturbation vector
		status = ranecu(seed);
		normalpert[2] = 2.0f*status - 1.0f;	// random number between (-1,1)
		rr = sqrt(1.0f - status*status);
		status = ranecu(seed);
		theta = status * 2.0f * pi;		// random number between (0,2pi)

		normalpert[0] = rr * cos(theta);
		normalpert[1] = rr * sin(theta);

		// normalize the perturbed vector
		normalize_base = sqrt( pow(normalpert[0],2) + pow(normalpert[1],2) + pow(normalpert[2],2) );
	
		normalpert[0] = normalpert[0]/normalize_base;
		normalpert[1] = normalpert[1]/normalize_base;
		normalpert[2] = normalpert[2]/normalize_base;

		// rough normal = beta*perturbed + original normal
		rough_normal[0] = beta * normalpert[0] + normal[0];	
		rough_normal[1] = beta * normalpert[1] + normal[1];
		rough_normal[2] = beta * normalpert[2] + normal[2];

		// normalize rough normal
		normalize_base = sqrt( pow(rough_normal[0],2) + pow(rough_normal[1],2) + pow(rough_normal[2],2) );

		normal[0] = rough_normal[0]/normalize_base; 
		normal[1] = rough_normal[1]/normalize_base;
		normal[2] = rough_normal[2]/normalize_base;

	return;
	}	

#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// determine if the photon gets detected
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
__device__ inline int detection(float3 *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, 
unsigned long long int *myimage, size_t pitch, struct start_info *info, int mytid, int *num_detected_primary, float sensorRefl, float d_min, int2* seed, 
float3 *dcos, float3 *normal, float bulk_abscoeff, float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, float *Xc, float *Yc,
struct histStruct * dPhotonHist, int boolToCollectData)
{
	int result = 0;		// flag - photon detected at sensor plane or absorbed in the bulk
	int absflag = 0;	// flag - photon absorbed in the bulk (yes=1; no=0); returned from isotropic
	int ii = 0, jj = 0;	// x,y pixel number for detected photons
		
	/*if(mytid == 6 || mytid == 7)
	{
			printf("tid: %d in detection\n");
	}*/
		
	// if z = -H/2 (sensor plane), photon detected
	if (fabs(pos->z - (float)(-H/2.0f)) < epsilon) 
	{
		if(ranecu(seed) < sensorRefl)	// non-ideal sensor - reflects back sensorRefl% of photons into the current column; detects rest
		{
			// normal pointing +z-direction
			normal->x = 0.0f; 
			normal->y = 0.0f; 
			normal->z = 1.0f;

			// reflected dir. cosine from the sensor plane (bottom surface) (bottom surface is smooth, no need to perturb the normal; specular reflection)
			trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, mytid, info); // reflection only, so 'refl_theta,trans_theta' = 0

			// using above calculated dir. cosine, move the photon within the column
			absflag = isotropic(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, &info[mytid], mynum_rebound, Xc, Yc, mytid, dPhotonHist, boolToCollectData);

			if(absflag == 1)	// absorbed in bulk
				result = 1;
			else
				result = 0;
		}
		else				// photon detected  
		{
			result = 1;
			atomicAdd(&num_detect, 1);	// increment photon detected counter

			if(mytid < dev_numPhotonHist && boolToCollectData == 1)
			{
				dPhotonHist[mytid].x[dPhotonHist[mytid].histCounter] = pos->x;
				dPhotonHist[mytid].y[dPhotonHist[mytid].histCounter] = pos->y;
				dPhotonHist[mytid].z[dPhotonHist[mytid].histCounter] = pos->z;
				dPhotonHist[mytid].Xc[dPhotonHist[mytid].histCounter] = *Xc;
				dPhotonHist[mytid].Yc[dPhotonHist[mytid].histCounter] = *Yc;
				dPhotonHist[mytid].terminated[dPhotonHist[mytid].histCounter] = 1;
				dPhotonHist[mytid].histCounter ++;
			}
			
			/*if(mytid == 6 || mytid == 7)
			{
				printf("detection - tid: %d, %f %f %f %f %f\n", mytid, pos->x, pos->y, pos->z, *Xc, *Yc);
			}*/
			
			ii = floor((pos->x-lbound_x)/pixelsize);	// x,y pixel number in the sensor plane of the detector
			jj = floor((pos->y-lbound_y)/pixelsize);

			// if the photon gets detected within lower and upper bounds of PRF: accumulate the signal contribution
			if( (pos->x <= ubound_x) && (pos->y <= ubound_y) && (pos->x >= lbound_x) && (pos->y >= lbound_y) )
			{	
				unsigned long long int* current_img = (unsigned long long int*)((char*)myimage + ii * pitch);
				atomicAdd(&current_img[jj],1);
			}
			atomicAdd(&num_detected_primary[info[mytid].str_histnum-1],1);	// increment # detected per primary
		} 
	}
	else
		result = 0;

	return result;
}	
#else	// C code

	int detection(float *pos, float H, int pixelsize, float lbound_x, float lbound_y, float ubound_x, float ubound_y, unsigned long long int *myimage, struct start_info info, 
	float sensorRefl, float d_min, int* seed, float *dcos, float *normal, float bulk_abscoeff, float R, float xdetector, float ydetector, unsigned long long int mynum_rebound, 
	int ydim, int *h_num_detected_prim)
	{

		int result = 0;		// flag - photon detected at sensor plane or absorbed in the bulk
		int absflag = 0;	// flag - photon absorbed in the bulk (yes=1; no=0); returned from isotropic
		int ii = 0, jj = 0;	// x,y pixel number for detected photons


		// if z = -H/2 (sensor plane), photon detected
		if (fabs(pos[2] - (float)(-H/2.0f)) < epsilon) 
		 {
			if(ranecu(seed) < sensorRefl)		// non-ideal sensor - reflects back sensorRefl% of photons into the current column; detects rest
			{

				// normal pointing +z-direction
				normal[0] = 0.0f; 
				normal[1] = 0.0f; 
				normal[2] = 1.0f;

				// reflected dir. cosine from the sensor plane (bottom surface) (bottom surface is smooth, no need to perturb the normal; specular reflection)
				trans_dir_cos(dcos, normal, 0.0f, 0.0f, 0, info); // reflection only, so 'refl_theta, trans_theta' = 0

				// using above calculated dir. cosine, move the photon within the column
				absflag = isotropic(pos, dcos, seed, bulk_abscoeff, R, H, xdetector, ydetector, info, mynum_rebound);

				if(absflag == 1)	// absorbed in bulk
					result = 1;
				else
					result = 0;

			}
			else					// photon detected  
			{
				result = 1;
				num_detect++;
		
				ii = floor((pos[0]-lbound_x)/pixelsize);	// x,y pixel number in the sensor plane of the detector
				jj = floor((pos[1]-lbound_y)/pixelsize);

				// if the photon gets detected within lower and upper bounds of PRF: accumulate the signal contribution
				if( (pos[0] <= ubound_x) && (pos[1] <= ubound_y) && (pos[0] >= lbound_x) && (pos[1] >= lbound_y) )
				 {	
					outputimage_.newimageopt[ii][jj]++;
				 }

				h_num_detected_prim[info.str_histnum]++;	// increment # detected per primary
			}
			 
		 }
		else
		    	result = 0;
	  
	 return result;
	}

#endif


////////////////////////////////////////////////////////////////////////////////
//! Initialize the pseudo-random number generator (PRNG) RANECU to a position
//! far away from the previous history (leap frog technique).
//!
//! Each calculated seed initiates a consecutive and disjoint sequence of
//! pseudo-random numbers with length LEAP_DISTANCE, that can be used to
//! in a parallel simulation (Sequence Splitting parallelization method).
//! The basic equation behind the algorithm is:
//!    S(i+j) = (a**j * S(i)) MOD m = [(a**j MOD m)*S(i)] MOD m  ,
//! which is described in:
//!   P L'Ecuyer, Commun. ACM 31 (1988) p.742
//!
//! This function has been adapted from "seedsMLCG.f", see:
//!   A Badal and J Sempau, Computer Physics Communications 175 (2006) p. 440-450
//!
//!       @param[in] history   Particle bach number.
//!       @param[in] seed_input   Initial PRNG seed input (used to initiate both MLCGs in RANECU).
//!       @param[out] seed   Initial PRNG seeds for the present history.
//!
////////////////////////////////////////////////////////////////////////////////
// -- Upper limit of the number of random values sampled in a single track:
#define  LEAP_DISTANCE    1000
// -- Multipliers and moduli for the two MLCG in RANECU:
#define  a1_RANECU       40014
#define  m1_RANECU  2147483563
#define  a2_RANECU       40692
#define  m2_RANECU  2147483399

#ifdef USING_CUDA
	__device__ inline void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int2* seed)
	{
	  // -- Move the RANECU generator to a unique position for the current batch of histories:
	  //    I have to use an "unsigned long long int" value to represent all the simulated histories in all previous batches
	  //    The maximum unsigned long long int value is ~1.8e19: if history >1.8e16 and LEAP_DISTANCE==1000, 'leap' will overflow.
	  // **** 1st MLCG:
	  unsigned long long int leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  int y = 1;
	  int z = a1_RANECU;
	  // -- Calculate the modulo power '(a^leap)MOD(m)' using a divide-and-conquer algorithm adapted to modulo arithmetic
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  Equivalent: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);  
	      y = abMODm(m1_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODm(m1_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm1 = y;                 // Exponentiation finished:  AjMODm = expMOD = y = a^j

	  // -- Compute and display the seeds S(i+j), from the present seed S(i), using the previously calculated value of (a^j)MOD(m):
	  //         S(i+j) = [(a**j MOD m)*S(i)] MOD m
	  //         S_i = abMODm(m,S_i,AjMODm)
	  seed->x = abMODm(m1_RANECU, seed_input, y);     // Using the input seed as the starting seed

	  // **** 2nd MLCG (repeating the previous calculation for the 2nd MLCG parameters):
	  leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  y = 1;
	  z = a2_RANECU;
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  Equivalent: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	      y = abMODm(m2_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODm(m2_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm2 = y;
	  seed->y = abMODm(m2_RANECU, seed_input, y);     // Using the input seed as the starting seed
	}
#else
	void init_PRNG(int history_batch, int histories_per_thread, int seed_input, int* seed)
	{
	  // -- Move the RANECU generator to a unique position for the current batch of histories:
	  //    I have to use an "unsigned long long int" value to represent all the simulated histories in all previous batches
	  //    The maximum unsigned long long int value is ~1.8e19: if history >1.8e16 and LEAP_DISTANCE==1000, 'leap' will overflow.
	  // **** 1st MLCG:
	  unsigned long long int leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  int y = 1;
	  int z = a1_RANECU;
	  // -- Calculate the modulo power '(a^leap)MOD(m)' using a divide-and-conquer algorithm adapted to modulo arithmetic
	  for(;;)
	  {
	      // printf(" leap, leap>>1, leap&1: %d, %d, %d\n",leap, leap>>1, leap&1);  

	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  !!DeBuG!! OLD: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2); 
	      y = abMODm(m1_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODm(m1_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm1 = y;                 // Exponentiation finished:  AjMODm = expMOD = y = a^j

	  // -- Compute and display the seeds S(i+j), from the present seed S(i), using the previously calculated value of (a^j)MOD(m):
	  //         S(i+j) = [(a**j MOD m)*S(i)] MOD m
	  //         S_i = abMODm(m,S_i,AjMODm)
	  seed[0] = abMODm(m1_RANECU, seed_input, y);     // Using the input seed as the starting seed

	  // **** 2nd MLCG (repeating the previous calculation for the 2nd MLCG parameters):
	  leap = ((unsigned long long int)(history_batch+1))*(histories_per_thread*LEAP_DISTANCE);
	  y = 1;
	  z = a2_RANECU;
	  for(;;)
	  {
	    // (A2) Halve n, and store the integer part and the residue
	    if (0!=(leap&01))  // (bit-wise operation for MOD(leap,2), or leap%2 ==> proceed if leap is an odd number)  !!DeBuG!! OLD: t=(short)(leap%2);
	    {
	      leap >>= 1;     // Halve n moving the bits 1 position right. Equivalent to:  leap=(leap/2); 
	      y = abMODm(m2_RANECU,z,y);      // (A3) Multiply y by z:  y = [z*y] MOD m
	      if (0==leap) break;         // (A4) leap==0? ==> finish
	    }
	    else           // (leap is even)
	    {
	      leap>>= 1;     // Halve leap moving the bits 1 position right. Equivalent to:  leap=(leap/2);
	    }
	    z = abMODm(m2_RANECU,z,z);        // (A5) Square z:  z = [z*z] MOD m
	  }
	  // AjMODm2 = y;
	  seed[1] = abMODm(m2_RANECU, seed_input, y);     // Using the input seed as the starting seed

	}
#endif


/////////////////////////////////////////////////////////////////////
//!  Calculate "(a1*a2) MOD m" with 32-bit integers and avoiding   **
//!  the possible overflow, using the Russian Peasant approach     **
//!  modulo m and the approximate factoring method, as described   **
//!  in:  L'Ecuyer and Cote, ACM Trans. Math. Soft. 17 (1991)      **
//!                                                                **
//!  This function has been adapted from "seedsMLCG.f", see:       **
//!  Badal and Sempau, Computer Physics Communications 175 (2006)  **
//!                                                                **
//!    Input:          0 < a1 < m                                  **
//!                    0 < a2 < m                                  **
//!                                                                **
//!    Return value:  (a1*a2) MOD m                                **
//!                                                                **
/////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline int abMODm(int m, int a, int s)
	{
	  // CAUTION: the input parameters are modified in the function but should not be returned to the calling function! (pass by value!)
	  int q, k;
	  int p = -m;            // p is always negative to avoid overflow when adding

	  // ** Apply the Russian peasant method until "a =< 32768":
	  while (a>32768)        // We assume '32' bit integers (4 bytes): 2^(('32'-2)/2) = 32768
	  {
	    if (0!=(a&1))        // Store 's' when 'a' is odd     Equivalent code:   if (1==(a%2))
	    {
	      p += s;
	      if (p>0) p -= m;
	    }
	    a >>= 1;             // Half a (move bits 1 position right)   Equivalent code: a = a/2;
	    s = (s-m) + s;       // Double s (MOD m)
	    if (s<0) s += m;     // (s is always positive)
	  }

	  // ** Employ the approximate factoring method (a is small enough to avoid overflow):
	  q = (int) m / a;
	  k = (int) s / q;
	  s = a*(s-k*q)-k*(m-q*a);
	  while (s<0)
	    s += m;

	  // ** Compute the final result:
	  p += s;
	  if (p<0) p += m;

	  return p;
	}
#else
	int abMODm(int m_par, int a_par, int s_par)
	{
	  // CAUTION: the input parameters are modified in the function but should not be returned to the calling function! (pass by value!)   !!DeBuG!!
	  int mval,aval,sval;
	  mval=m_par; aval=a_par; sval=s_par;
	  
	  int qval, kval;
	  int pval = -mval;            // p is always negative to avoid overflow when adding

	  // ** Apply the Russian peasant method until "a =< 32768":
	  while (aval>32768)        // We assume '32' bit integers (4 bytes): 2^(('32'-2)/2) = 32768
	  {
	    if (0!=(aval&1))        // Store 's' when 'a' is odd    !!DeBuG!! OLD code:   if (1==(a%2))
	    {
	      pval += sval;
	      if (pval>0) pval -= mval;
	    }
	    aval >>= 1;             // Half a (move bits 1 position right)        
	    sval = (sval-mval) + sval;       // float s (MOD m)
	    if (sval<0) sval += mval;     // (s is always positive)
	  }

	  // ** Employ the approximate factoring method (a is small enough to avoid overflow):
	  qval = (int) mval / aval;
	  kval = (int) sval / qval;
	  sval = aval*(sval-kval*qval)-kval*(mval-qval*aval);
	  while (sval<0)
	    sval += mval;

	  // ** Compute the final result:
	  pval += sval;
	  if (pval<0) pval += mval;

	  return pval;
	}
#endif


////////////////////////////////////////////////////////////////////////////////
//! Pseudo-random number generator (PRNG) RANECU returning a float value
//! (single precision version).
//!
//!       @param[in,out] seed   PRNG seed (seed kept in the calling function and updated here).
//!       @return   PRN double value in the open interval (0,1)
//!
////////////////////////////////////////////////////////////////////////////////

#ifdef USING_CUDA
	__device__ inline float ranecu(int2* seed)
	{
	//return (float(seed->x%100)*0.01f+0.005f)  ;

	  int i1 = (int)(seed->x/53668);
	  seed->x = 40014*(seed->x-i1*53668)-i1*12211;

	  int i2 = (int)(seed->y/52774);
	  seed->y = 40692*(seed->y-i2*52774)-i2*3791;

	  if (seed->x < 0) seed->x += 2147483563;
	  if (seed->y < 0) seed->y += 2147483399;

	  i2 = seed->x-seed->y;
	  if (i2 < 1) i2 += 2147483562;

	  return (__int2float_rn(i2)*4.65661305739e-10f);        // 4.65661305739e-10 == 1/2147483563

	}
#else
	float ranecu(int* seed)
	{
	  int i1 = (int)(seed[0]/53668);
	  seed[0] = 40014*(seed[0]-i1*53668)-i1*12211;

	  int i2 = (int)(seed[1]/52774);
	  seed[1] = 40692*(seed[1]-i2*52774)-i2*3791;

	  if (seed[0] < 0) seed[0] += 2147483563;
	  if (seed[1] < 0) seed[1] += 2147483399;

	  i2 = seed[0]-seed[1];
	  if (i2 < 1) i2 += 2147483562;

	  const float USCALE = 1.0/2147483563.0;       
	  return ((float)(i2*USCALE));

	}
#endif

