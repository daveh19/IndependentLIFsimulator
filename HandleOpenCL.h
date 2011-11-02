/*
 *  HandleOpenCL.h
 *  XclNet
 *
 *  Created by David Higgins on 19/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"


// Use a static data size for simplicity
//
#define DATA_SIZE (8)

// Wrap variables in a struct
typedef struct {
	int err;                            // error code returned from api calls
	
    //float data[DATA_SIZE];              // original data set given to device
    //float results[DATA_SIZE];           // results returned from device
	
    size_t global;                      // global domain size for our calculation: work size per dimension
    size_t local;                       // local domain size for our calculation: work-group size per dimension
	
    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
	
    //cl_mem input;                       // device memory object used for the input array
    //cl_mem output;                      // device memory object used for the output array
	
	cl_mem input_v;
	cl_mem input_current;
	cl_mem input_gauss;
	cl_mem input_spike;
	
	cl_mem output_v;
	cl_mem output_spike;
	
} CL;

char* readKernelSource(char * filename);
int connectToComputeDevice(CL *cl);
int createComputeContext(CL *cl);
int createCommandQueue(CL *cl);
int createProgram(CL *cl, char ** KernelSource);
int buildProgram(CL *cl);
int createKernel(CL *cl, char * k_name);
//int createIObufs(CL *cl, unsigned int count);
int createLifIObufs(CL *cl);
//int enqueueInputBuf(CL *cl, unsigned int count);
int enqueueLifInputBuf(CL *cl, cl_LIFNeuron *lif);
//int setKernelArgs(CL *cl, unsigned int count);
int setLifKernelArgs(CL *cl, cl_LIFNeuron *lif);
int getMaxWorkSize(CL *cl);
//int enqueueKernel(CL *cl, unsigned int count);
int enqueueLifKernel(CL *cl);
void waitForKernel(CL *cl);
//int enqueueOutputBuf(CL *cl, unsigned int count);
int enqueueLifOutputBuf(CL *cl, cl_LIFNeuron *lif);

//void shutdownKernel(CL *cl);
void shutdownLifKernel(CL *cl);

// My Macro functions
int setupCL(CL *cl);
int makeProgram(CL *cl, char* KernelSource, char* k_name);
