/*
 *  HandleOpenCL.h
 *  XclNet
 *
 *  Created by David Higgins on 19/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralIncludes.h"


// Use a static data size for simplicity
//
#define DATA_SIZE (1024)

// Wrap variables in a struct
typedef struct {
	int err;                            // error code returned from api calls
	
    float data[DATA_SIZE];              // original data set given to device
    float results[DATA_SIZE];           // results returned from device
	
    size_t global;                      // global domain size for our calculation: work size per dimension
    size_t local;                       // local domain size for our calculation: work-group size per dimension
	
    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
	
    cl_mem input;                       // device memory object used for the input array
    cl_mem output;                      // device memory object used for the output array
	
	
} CL;

char* readKernelSource(char * filename);
int connectToComputeDevice(CL *cl);
int createComputeContext(CL *cl);
int createCommandQueue(CL *cl);
int createProgram(CL *cl, char ** KernelSource);
int buildProgram(CL *cl);
int createKernel(CL *cl, char * k_name);
int createIObufs(CL *cl, unsigned int count);
int enqueueInputBuf(CL *cl, unsigned int count);
int setKernelArgs(CL *cl, unsigned int count);
int getMaxWorkSize(CL *cl);
int enqueueKernel(CL *cl, unsigned int count);
void waitForKernel(CL *cl);
int enqueueOutputBuf(CL *cl, unsigned int count);

void shutdownKernel(CL *cl);

// My Macro functions
int setupCL(CL *cl);
int makeProgram(CL *cl, char* KernelSource, char* k_name);
