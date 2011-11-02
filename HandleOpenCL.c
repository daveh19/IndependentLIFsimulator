/*
 *  HandleOpenCL.c
 *  XclNet
 *
 *  Created by David Higgins on 19/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "HandleOpenCL.h"


char* readKernelSource(char * filename){
	// Try loading the kernel from a source file
	FILE *f_kernel;
	
	f_kernel = fopen(filename, "rb");
	if (f_kernel == NULL){
		perror("Error in loading kernel from file");
	}
	
	fseek(f_kernel, 0, SEEK_END);
	long pos = ftell(f_kernel);
	fseek(f_kernel, 0, SEEK_SET);
	
	char *KernelSource = malloc(pos);
	fread(KernelSource, pos, 1, f_kernel);
	fclose(f_kernel);
	printf("Source:\n %s", KernelSource);	
	
	return KernelSource;
}

int setupCL(CL *cl){
	// Create initial OpenCL context and queue
	//
	
	if( connectToComputeDevice(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( createComputeContext(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( createCommandQueue(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	return !(EXIT_FAILURE);
}

int connectToComputeDevice(CL *cl){
	// Connect to a compute device
	//
	int gpu = 1;
	
	printf("connecting to compute device...\n");
	
	(*cl).err = clGetDeviceIDs(NULL, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &(*cl).device_id, NULL);
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to create a device group!\n");
		return EXIT_FAILURE;
	}
	//printf("test %d\n", !(EXIT_FAILURE));
	return !(EXIT_FAILURE);
}

int createComputeContext(CL *cl){
	// Create a compute context
    //
	
	printf("creating compute context...\n");
	
    (*cl).context = clCreateContext(0, 1, &(*cl).device_id, NULL, NULL, &(*cl).err);
    if (!(*cl).context)
    {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int createCommandQueue(CL *cl){
	// Create a command commands
    //
	
	printf("creating command queue...\n");

    (*cl).commands = clCreateCommandQueue((*cl).context, (*cl).device_id, 0, &(*cl).err);
    if (!(*cl).commands)
    {
        printf("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int makeProgram(CL *cl, char* KernelSource, char* k_name){
	// Make a compute kernel from source
	//
	
	//printf("Making program from kernel source:\n");
	//printf("%s", KernelSource);
	
	if( createProgram(cl, &KernelSource) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( buildProgram(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( createKernel(cl, k_name) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	return !(EXIT_FAILURE);
}

int createProgram(CL *cl, char ** KernelSource){
    // Create the compute program from the source buffer
    //
	
	printf("creating program from source...\n");
	
	printf("%s", *KernelSource);
	
	//printf("Kernel source: %s\n", *KernelSource);
	
    (*cl).program = clCreateProgramWithSource((*cl).context, 1, (const char **) KernelSource, NULL, &(*cl).err);
    if (!(*cl).program)
    {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int buildProgram(CL *cl){
	// Build the program executable
	//
	
	printf("building program...\n");
	
	(*cl).err = clBuildProgram((*cl).program, 0, NULL, NULL, NULL, NULL);
	if ((*cl).err != CL_SUCCESS)
	{
		size_t len;
		char buffer[2048];
		
		printf("Error: Failed to build program executable!\n");
		clGetProgramBuildInfo((*cl).program, (*cl).device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

int createKernel(CL *cl, char * k_name){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("creating kernel (%s)...\n", k_name);
	
	(*cl).kernel = clCreateKernel((*cl).program, k_name, &(*cl).err);
	if (!(*cl).kernel || (*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to create compute kernel!\n");
		exit(1);
	}
	return !(EXIT_FAILURE);
}

int createIObufs(CL *cl, unsigned int count){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("creating io buffers...\n");
	
    // Create the input and output arrays in device memory for our calculation
    //
    (*cl).input = clCreateBuffer((*cl).context,  CL_MEM_READ_ONLY,  sizeof(float) * count, NULL, NULL);
    (*cl).output = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(float) * count, NULL, NULL);
    if (!(*cl).input || !(*cl).output)
    {
        printf("Error: Failed to allocate device memory!\n");
	    exit(1);
	}
	return !(EXIT_FAILURE);
}

int enqueueInputBuf(CL *cl, unsigned int count){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("enqueueing input buffer...\n");
	
    // Write our data set into the input array in device memory
    //
    (*cl).err = clEnqueueWriteBuffer((*cl).commands, (*cl).input, CL_TRUE, 0, sizeof(float) * count, (*cl).data, 0, NULL, NULL);
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
	return !(EXIT_FAILURE);
}

int setKernelArgs(CL *cl, unsigned int count){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("setting args for compute kernel...\n");
	
    // Set the arguments to our compute kernel
    //
    (*cl).err = 0;
    (*cl).err  = clSetKernelArg((*cl).kernel, 0, sizeof(cl_mem), &(*cl).input);
    (*cl).err |= clSetKernelArg((*cl).kernel, 1, sizeof(cl_mem), &(*cl).output);
    (*cl).err |= clSetKernelArg((*cl).kernel, 2, sizeof(unsigned int), &count);
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %d\n", (*cl).err);
        exit(1);
    }
	return !(EXIT_FAILURE);
}

int getMaxWorkSize(CL *cl){
	// Get the maximum work group size for executing the kernel on the device
	//
	
	printf("getting max work group size...\n");
	
	(*cl).err = clGetKernelWorkGroupInfo((*cl).kernel, (*cl).device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof((*cl).local), &(*cl).local, NULL);
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to retrieve kernel work group info! %d\n", (*cl).err);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

int enqueueKernel(CL *cl, unsigned int count){
	// Execute the kernel over the entire range of our 1d input data set
	// using the maximum number of work group items for this device
	//
	
	printf("sending the kernel to the process queue...\n");
	
	(*cl).global = count;
	(*cl).err = clEnqueueNDRangeKernel((*cl).commands, (*cl).kernel, 1, NULL, &(*cl).global, &(*cl).local, 0, NULL, NULL);
	if ((*cl).err)
	{
		printf("Error: Failed to execute kernel!\n");
		return EXIT_FAILURE;
	}
	return !(EXIT_FAILURE);
}

void waitForKernel(CL *cl){
	// Wait for the command commands to get serviced before reading back results
	//
	
	printf("waiting for kernel to finish...\n");
	
	clFinish((*cl).commands);
}

int enqueueOutputBuf(CL *cl, unsigned int count){
	// Read back the results from the device to verify the output
	//
	
	printf("reading output from kernel...\n");
	
	(*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).output, CL_TRUE, 0, sizeof(float) * count, (*cl).results, 0, NULL, NULL );
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to read output array! %d\n", (*cl).err);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

void shutdownKernel(CL *cl){
	// Shutdown and cleanup
	//
	
	printf("shutting down and cleaning up kernel memory...");
	
	clReleaseMemObject((*cl).input);
	clReleaseMemObject((*cl).output);
	clReleaseProgram((*cl).program);
	clReleaseKernel((*cl).kernel);
	clReleaseCommandQueue((*cl).commands);
	clReleaseContext((*cl).context);
	
	printf("done\n");
}