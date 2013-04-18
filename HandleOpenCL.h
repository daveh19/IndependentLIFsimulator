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
#include "cl_Synapse.h"

//int job_size;

typedef struct random_struct_marsaglia{
	//float value;
	unsigned int * d_z;
	unsigned int * d_w;
	unsigned int * d_jsr;
	unsigned int * d_jcong;
} cl_MarsagliaStruct;

// Wrap variables in a struct
typedef struct cl_struct{
	int err;                            // error code returned from api calls
	
    size_t global;                      // global domain size for our calculation: work size per dimension
    size_t local;                       // local domain size for our calculation: work-group size per dimension
	
	cl_platform_id platform;	
    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
	
	
	cl_mem gauss;
	
	// LIF specific memory streams
	cl_mem input_v;
	cl_mem input_current;
	cl_mem input_spike;
	
	// Synapse specific memory streams
	cl_mem rho;
	cl_mem ca;
	cl_mem pre_spike;
	cl_mem post_spike;
	
	
	// Random number memory streams
	cl_mem d_z;
	cl_mem d_w;
	cl_mem d_jsr;
	cl_mem d_jcong;
	
	
	int job_size;
} CL;


//static char *print_cl_errstring(cl_int err);
char* readKernelSource(char * filename);
int getPlatformIDs(CL *cl);
int connectToComputeDevice(CL *cl);
int createComputeContext(CL *cl);
int createCommandQueue(CL *cl);
int createProgram(CL *cl, char ** KernelSource);
int buildProgram(CL *cl);
int createKernel(CL *cl, char * k_name);
//int createIObufs(CL *cl, unsigned int count);
int createLifIObufs(CL *cl);
int createSynIObufs(CL *cl);
//int enqueueInputBuf(CL *cl, unsigned int count);
int enqueueLifInputBuf(CL *cl, cl_LIFNeuron *lif, cl_MarsagliaStruct *rnd);
int enqueueSynInputBuf(CL *cl, cl_Synapse *syn, SynapseConsts *syn_const, cl_MarsagliaStruct *rnd);
//int setKernelArgs(CL *cl, unsigned int count);
int setLifKernelArgs(CL *cl, cl_LIFNeuron *lif);
int setSynKernelArgs(CL *cl, cl_Synapse *syn, SynapseConsts *syn_const);
int getMaxWorkSize(CL *cl);
//int enqueueKernel(CL *cl, unsigned int count);
int enqueueLifKernel(CL *cl);
int enqueueSynKernel(CL *cl);
void waitForKernel(CL *cl);
//int enqueueOutputBuf(CL *cl, unsigned int count);
int enqueueLifOutputBuf(CL *cl, cl_LIFNeuron *lif, cl_MarsagliaStruct *rnd);
int enqueueSynOutputBuf(CL *cl, cl_Synapse *syn, SynapseConsts *syn_const, cl_MarsagliaStruct *rnd);

//void shutdownKernel(CL *cl);
void shutdownLifKernel(CL *cl);
void shutdownSynKernel(CL *cl);

// My Macro functions
int setupCL(CL *cl);
int makeProgram(CL *cl, char* KernelSource, char* k_name);
