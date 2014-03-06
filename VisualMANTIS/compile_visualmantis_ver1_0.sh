#!/bin/bash

# CUDA codes
nvcc -O3 -use_fast_math -DUSING_CUDA -c hybridMANTIS_cuda_ver1_0_LB.cu -I/usr/local/cuda/include -I/home/han/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/local/cuda/lib64/ -lcutil_x86_64 -lcudart -lgsl -lgslcblas -lm -arch sm_20 --ptxas-options=-v

nvcc -O3 -use_fast_math -DUSING_CUDA -c hybridMANTIS_cuda_ver1_0.cu -I/usr/local/cuda/include -I/home/han/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/local/cuda/lib64/ -lcutil_x86_64 -lcudart -lgsl -lgslcblas -lm -arch sm_20 --ptxas-options=-v

nvcc -O3 -use_fast_math -DUSING_CUDA -c visualMANTIS_cuda_ver1_0.cu -I/usr/local/cuda/include -I/home/han/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/local/cuda/lib64/ -lcutil_x86_64 -lcudart -lgsl -lgslcblas -lm -arch sm_20 --ptxas-options=-v

# C codes
gcc -I/usr/include -c hybridMANTIS_c_ver1_0_LB.c -lgsl -lgslcblas -lm

gcc -I/usr/include -c hybridMANTIS_c_ver1_0.c -lgsl -lgslcblas -lm

# link with PENELOPE
gfortran -O3 -g hybridMANTIS_penEasy.f hybridMANTIS_tallyEnergyDepositionEvents.f penelope.f pengeom.f penvared.f hybridMANTIS_cuda_ver1_0_LB.o hybridMANTIS_cuda_ver1_0.o hybridMANTIS_c_ver1_0_LB.o hybridMANTIS_c_ver1_0.o -o hybridMANTIS_ver1_0.x -I/usr/local/cuda/include -I/home/han/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/lib -L/home/han/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/local/cuda/lib64 -lstdc++ -lcutil_x86_64 -lcudart -lgsl -lgslcblas -lm

# link with PENELOPE
gfortran -O3 -g hybridMANTIS_penEasy.f hybridMANTIS_tallyEnergyDepositionEvents.f penelope.f pengeom.f penvared.f hybridMANTIS_cuda_ver1_0_LB.o visualMANTIS_cuda_ver1_0.o hybridMANTIS_c_ver1_0_LB.o hybridMANTIS_c_ver1_0.o -o visualMANTIS_ver1_0.x -I/usr/local/cuda/include -I/home/han/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/lib -L/home/han/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/local/cuda/lib64 -lstdc++ -lcutil_x86_64 -lcudart -lgsl -lgslcblas -lm

#
fltk-config --use-gl --use-images --use-glut --compile visualization/visualizationGUI.cxx
#g++ -I/usr/local/include -I/usr/include/freetype2 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_THREAD_SAFE -D_REENTRANT -o 'visualizationGUI' 'visualizationGUI.cxx' /usr/local/lib/libfltk_gl.a -lGLU -lGL /usr/local/lib/libfltk_images.a -lpng -lz -ljpeg /usr/local/lib/libfltk.a -lXext -lXft -lfontconfig -lXinerama -lpthread -ldl -lm -lX11

cp visualizationGUI visualMANTIS_ver1_0.x example/

