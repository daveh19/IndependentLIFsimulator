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
	int gpu = USE_GPU;
	
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

//int createIObufs(CL *cl, unsigned int count){
//	// Create the compute kernel in the program we wish to run
//	//
//	
//	printf("creating io buffers...\n");
//	
//    // Create the input and output arrays in device memory for our calculation
//    //
//    (*cl).input = clCreateBuffer((*cl).context,  CL_MEM_READ_ONLY,  sizeof(float) * count, NULL, NULL);
//    (*cl).output = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(float) * count, NULL, NULL);
//    if (!(*cl).input || !(*cl).output)
//    {
//        printf("Error: Failed to allocate device memory!\n");
//	    exit(1);
//	}
//	return !(EXIT_FAILURE);
//}

int createLifIObufs(CL *cl){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("creating LIF io buffers...\n");
	
    // Create the input and output arrays in device memory for our calculation
    //
    (*cl).input_v = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_LIFS, NULL, NULL);
    (*cl).input_current = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_LIFS, NULL, NULL);
	(*cl).input_gauss = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_LIFS, NULL, NULL);
	(*cl).input_spike = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(unsigned int) * NO_LIFS, NULL, NULL);
	
	//(*cl).output_v = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(float) * NO_LIFS, NULL, NULL);
	//(*cl).output_spike = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(unsigned int) * NO_LIFS, NULL, NULL);
    if (!(*cl).input_v || !(*cl).input_current || !(*cl).input_gauss || !(*cl).input_spike /*|| !(*cl).output_v || !(*cl).output_spike*/)
    {
        printf("Error: Failed to allocate device memory!\n");
	    exit(1);
	}
	return !(EXIT_FAILURE);
}

int createSynIObufs(CL *cl){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("creating Synapse io buffers...\n");
	
    // Create the input and output arrays in device memory for our calculation
    //
    (*cl).rho = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_SYNS, NULL, NULL);
    (*cl).ca = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_SYNS, NULL, NULL);
	(*cl).input_gauss = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(float) * NO_SYNS, NULL, NULL);
	(*cl).pre_spike = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(unsigned int) * NO_SYNS, NULL, NULL);
	(*cl).post_spike = clCreateBuffer((*cl).context,  CL_MEM_READ_WRITE,  sizeof(unsigned int) * NO_SYNS, NULL, NULL);
	
	//(*cl).output_rho = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(float) * NO_SYNS, NULL, NULL);
	//(*cl).output_ca = clCreateBuffer((*cl).context, CL_MEM_WRITE_ONLY, sizeof(float) * NO_SYNS, NULL, NULL);
    if (!(*cl).rho || !(*cl).ca || !(*cl).input_gauss || !(*cl).pre_spike || !(*cl).post_spike /*|| !(*cl).output_rho || !(*cl).output_ca*/)
    {
        printf("Error: Failed to allocate device memory!\n");
	    exit(1);
	}
	return !(EXIT_FAILURE);
}

//int enqueueInputBuf(CL *cl, unsigned int count){
//	// Create the compute kernel in the program we wish to run
//	//
//	
//	printf("enqueueing input buffer...\n");
//	
//    // Write our data set into the input array in device memory
//    //
//	(*cl).err = clEnqueueWriteBuffer((*cl).commands, (*cl).input, CL_TRUE, 0, sizeof(float) * count, (*cl).data, 0, NULL, NULL);
//	if ((*cl).err != CL_SUCCESS)
//    {
//        printf("Error: Failed to write to source array!\n");
//        exit(1);
//    }
//	return !(EXIT_FAILURE);
//}

int enqueueLifInputBuf(CL *cl, cl_LIFNeuron *lif){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("enqueueing LIF input buffer...\n");
	
    // Write our data set into the input array in device memory
    //
    (*cl).err = clEnqueueWriteBuffer((*cl).commands, (*cl).input_v, CL_TRUE, 0, sizeof(float) * NO_LIFS, (*lif).V, 0, NULL, NULL);
    (*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).input_current, CL_TRUE, 0, sizeof(float) * NO_LIFS, (*lif).I, 0, NULL, NULL);
	(*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).input_gauss, CL_TRUE, 0, sizeof(float) * NO_LIFS, (*lif).gauss, 0, NULL, NULL);
	(*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).input_spike, CL_TRUE, 0, sizeof(unsigned int) * NO_LIFS, (*lif).time_since_spike, 0, NULL, NULL);
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
	return !(EXIT_FAILURE);
}

int enqueueSynInputBuf(CL *cl, cl_Synapse *syn){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("enqueueing Synapse input buffer...\n");
	
    // Write our data set into the input array in device memory
    //
    (*cl).err = clEnqueueWriteBuffer((*cl).commands, (*cl).rho, CL_TRUE, 0, sizeof(float) * NO_SYNS, (*syn).rho, 0, NULL, NULL);
    (*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).ca, CL_TRUE, 0, sizeof(float) * NO_SYNS, (*syn).ca, 0, NULL, NULL);
	(*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).input_gauss, CL_TRUE, 0, sizeof(float) * NO_SYNS, (*syn).gauss, 0, NULL, NULL);
	(*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).pre_spike, CL_TRUE, 0, sizeof(unsigned int) * NO_SYNS, (*syn).preT[0], 0, NULL, NULL);
    (*cl).err |= clEnqueueWriteBuffer((*cl).commands, (*cl).post_spike, CL_TRUE, 0, sizeof(unsigned int) * NO_SYNS, (*syn).postT, 0, NULL, NULL);
	if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
	return !(EXIT_FAILURE);
}

//int setKernelArgs(CL *cl, unsigned int count){
//	// Create the compute kernel in the program we wish to run
//	//
//	
//	printf("setting args for compute kernel...\n");
//	
//    // Set the arguments to our compute kernel
//    //
//    (*cl).err = 0;
//    (*cl).err  = clSetKernelArg((*cl).kernel, 0, sizeof(cl_mem), &(*cl).input);
//    (*cl).err |= clSetKernelArg((*cl).kernel, 1, sizeof(cl_mem), &(*cl).output);
//    (*cl).err |= clSetKernelArg((*cl).kernel, 2, sizeof(unsigned int), &count);
//    if ((*cl).err != CL_SUCCESS)
//    {
//        printf("Error: Failed to set kernel arguments! %d\n", (*cl).err);
//        exit(1);
//    }
//	return !(EXIT_FAILURE);
//}

int setLifKernelArgs(CL *cl, cl_LIFNeuron *lif){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("setting args for LIF compute kernel...\n");
	
    // Set the arguments to our compute kernel
    //
    (*cl).err = 0;
    (*cl).err  = clSetKernelArg((*cl).kernel, 0, sizeof(cl_mem), &(*cl).input_v);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 1, sizeof(cl_mem), &(*cl).input_current);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 2, sizeof(cl_mem), &(*cl).input_gauss);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 3, sizeof(cl_mem), &(*cl).input_spike);
	
	//(*cl).err  |= clSetKernelArg((*cl).kernel, 4, sizeof(cl_mem), &(*cl).output_v);
	//(*cl).err  |= clSetKernelArg((*cl).kernel, 5, sizeof(cl_mem), &(*cl).output_spike);
	
	(*cl).err  |= clSetKernelArg((*cl).kernel, 4, sizeof(float), &(*lif).v_rest);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 5, sizeof(float), &(*lif).v_reset);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 6, sizeof(float), &(*lif).v_threshold);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 7, sizeof(float), &(*lif).r_m);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 8, sizeof(float), &(*lif).c_m);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 9, sizeof(float), &(*lif).sigma);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 10, sizeof(float), &(*lif).refrac_time);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 11, sizeof(float), &(*lif).dt);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 12, sizeof(unsigned int), &(*lif).no_lifs);
	
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %d\n", (*cl).err);
        exit(1);
    }
	return !(EXIT_FAILURE);
}

int setSynKernelArgs(CL *cl, cl_Synapse *syn){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("setting args for Synapse compute kernel...\n");
	
    // Set the arguments to our compute kernel
    //
    (*cl).err = 0;
    (*cl).err  = clSetKernelArg((*cl).kernel, 0, sizeof(cl_mem), &(*cl).rho);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 1, sizeof(cl_mem), &(*cl).ca);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 2, sizeof(cl_mem), &(*cl).input_gauss);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 3, sizeof(cl_mem), &(*cl).pre_spike);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 4, sizeof(cl_mem), &(*cl).post_spike);
	
	//(*cl).err  |= clSetKernelArg((*cl).kernel, 5, sizeof(cl_mem), &(*cl).output_rho);
	//(*cl).err  |= clSetKernelArg((*cl).kernel, 6, sizeof(cl_mem), &(*cl).output_ca);
	
	(*cl).err  |= clSetKernelArg((*cl).kernel, 5, sizeof(float), &(*syn).gamma_p);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 6, sizeof(float), &(*syn).gamma_d);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 7, sizeof(float), &(*syn).theta_p);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 8, sizeof(float), &(*syn).theta_p);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 9, sizeof(float), &(*syn).tau);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 10, sizeof(float), &(*syn).tau_ca);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 11, sizeof(float), &(*syn).c_pre);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 12, sizeof(float), &(*syn).c_post);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 13, sizeof(float), &(*syn).sigma);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 14, sizeof(float), &(*syn).dt);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 15, sizeof(unsigned int), &(*syn).no_syns);
	
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
	//size_t my_local_var;
	
	printf("getting max work group size...\n");
	
	(*cl).err = clGetKernelWorkGroupInfo((*cl).kernel, (*cl).device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof((*cl).local), &(*cl).local, NULL);
	//(*cl).err = clGetKernelWorkGroupInfo((*cl).kernel, (*cl).device_id, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(my_local_var), &my_local_var, NULL);//Doesn't work under apple
	//(*cl).err = clGetKernelWorkGroupInfo((*cl).kernel, (*cl).device_id, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(my_local_var), &my_local_var, NULL);
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to retrieve kernel work group info! %d\n", (*cl).err);
		exit(1);
	}
	printf("CL_KERNEL_WORK_GROUP_SIZE: %d\n", (int)(*cl).local);
	//printf("CL_KERNEL_LOCAL_MEM_SIZE: %ld\n", (long)my_local_var);
	return !(EXIT_FAILURE);
}

//int enqueueKernel(CL *cl, unsigned int count){
//	// Execute the kernel over the entire range of our 1d input data set
//	// using the maximum number of work group items for this device
//	//
//	
//	printf("sending the kernel to the process queue...\n");
//	
//	(*cl).global = count;
//	(*cl).err = clEnqueueNDRangeKernel((*cl).commands, (*cl).kernel, 1, NULL, &(*cl).global, &(*cl).local, 0, NULL, NULL);
//	if ((*cl).err)
//	{
//		printf("Error: Failed to execute kernel!\n");
//		return EXIT_FAILURE;
//	}
//	return !(EXIT_FAILURE);
//}

int enqueueLifKernel(CL *cl){
	// Execute the kernel over the entire range of our 1d input data set
	// using the maximum number of work group items for this device
	//
	
	printf("sending the LIF kernel to the process queue...\n");
	if((*cl).err){
		printf("Error already occurred\n");
	}
	(*cl).global = NO_LIFS;
	(*cl).local = fmin((*cl).local, (*cl).global); // Copes with case global < local
	while( (*cl).global % (*cl).local != 0){
		// Pad the global number of work items such that it is divided evenly by local
		(*cl).global++;
		//printf("New value for global: %d\n", (*cl).global);
	}
	printf("Executing with global: %d, local: %d, real no jobs: %d\n", (int)(*cl).global, (int)(*cl).local, NO_LIFS);

	(*cl).err = clEnqueueNDRangeKernel((*cl).commands, (*cl).kernel, 1, NULL, &(*cl).global, &(*cl).local, 0, NULL, NULL);
	if ((*cl).err)
	{
		printf("Error: Failed to execute kernel!\n");
		return EXIT_FAILURE;
	}
	return !(EXIT_FAILURE);
}

int enqueueSynKernel(CL *cl){
	// Execute the kernel over the entire range of our 1d input data set
	// using the maximum number of work group items for this device
	//
	
	printf("sending the Synapse kernel to the process queue...\n");
	if((*cl).err){
		printf("Error already occurred\n");
	}
	(*cl).global = NO_SYNS;
	(*cl).local = fmin((*cl).local, (*cl).global); // Copes with case global < local
	while( (*cl).global % (*cl).local != 0){
		// Pad the global number of work items such that it is divided evenly by local
		(*cl).global++;
		//printf("New value for global: %d\n", (*cl).global);
	}
	printf("Executing with global: %d, local: %d, real no jobs: %d\n", (int)(*cl).global, (int)(*cl).local, NO_SYNS);

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

//int enqueueOutputBuf(CL *cl, unsigned int count){
//	// Read back the results from the device to verify the output
//	//
//	
//	printf("reading output from kernel...\n");
//	
//	(*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).output, CL_TRUE, 0, sizeof(float) * count, (*cl).results, 0, NULL, NULL );
//	if ((*cl).err != CL_SUCCESS)
//	{
//		printf("Error: Failed to read output array! %d\n", (*cl).err);
//		exit(1);
//	}
//	return !(EXIT_FAILURE);
//}

int enqueueLifOutputBuf(CL *cl, cl_LIFNeuron *lif){
	// Read back the results from the device to verify the output
	//
	
	printf("reading output from LIF kernel...\n");
	
	//(*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).output_v, CL_TRUE, 0, sizeof(float) * NO_LIFS, (*lif).V, 0, NULL, NULL );
	//(*cl).err |= clEnqueueReadBuffer( (*cl).commands, (*cl).output_spike, CL_TRUE, 0, sizeof(unsigned int) * NO_LIFS, (*lif).time_since_spike, 0, NULL, NULL );
	(*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).input_v, CL_TRUE, 0, sizeof(float) * NO_LIFS, (*lif).V, 0, NULL, NULL );
	(*cl).err |= clEnqueueReadBuffer( (*cl).commands, (*cl).input_spike, CL_TRUE, 0, sizeof(unsigned int) * NO_LIFS, (*lif).time_since_spike, 0, NULL, NULL );
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to read output array! %d\n", (*cl).err);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

int enqueueSynOutputBuf(CL *cl, cl_Synapse *syn){
	// Read back the results from the device to verify the output
	//
	
	printf("reading output from Synapse kernel...\n");
	
	(*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).rho, CL_TRUE, 0, sizeof(float) * NO_SYNS, (*syn).rho, 0, NULL, NULL );
	(*cl).err |= clEnqueueReadBuffer( (*cl).commands, (*cl).ca, CL_TRUE, 0, sizeof(float) * NO_SYNS, (*syn).ca, 0, NULL, NULL );
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to read output array! %d\n", (*cl).err);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

//void shutdownKernel(CL *cl){
//	// Shutdown and cleanup
//	//
//	
//	printf("shutting down and cleaning up kernel memory...");
//	
//	clReleaseMemObject((*cl).input);
//	clReleaseMemObject((*cl).output);
//	clReleaseProgram((*cl).program);
//	clReleaseKernel((*cl).kernel);
//	clReleaseCommandQueue((*cl).commands);
//	clReleaseContext((*cl).context);
//	
//	printf("done\n");
//}

void shutdownLifKernel(CL *cl){
	// Shutdown and cleanup
	//
	
	printf("shutting down and cleaning up LIF kernel memory...");
	
	clReleaseMemObject((*cl).input_v);
	clReleaseMemObject((*cl).input_current);
	clReleaseMemObject((*cl).input_gauss);
	clReleaseMemObject((*cl).input_spike);
	//clReleaseMemObject((*cl).output_v);
	//clReleaseMemObject((*cl).output_spike);
	clReleaseProgram((*cl).program);
	clReleaseKernel((*cl).kernel);
	clReleaseCommandQueue((*cl).commands);
	clReleaseContext((*cl).context);
	
	printf("done\n");
}

void shutdownSynKernel(CL *cl){
	// Shutdown and cleanup
	//
	
	printf("shutting down and cleaning up Synapse kernel memory...");
	
	clReleaseMemObject((*cl).rho);
	clReleaseMemObject((*cl).ca);
	clReleaseMemObject((*cl).input_gauss);
	clReleaseMemObject((*cl).pre_spike);
	clReleaseMemObject((*cl).post_spike);
	//clReleaseMemObject((*cl).output_rho);
	//clReleaseMemObject((*cl).output_ca);
	clReleaseProgram((*cl).program);
	clReleaseKernel((*cl).kernel);
	clReleaseCommandQueue((*cl).commands);
	clReleaseContext((*cl).context);
	
	printf("done\n");
}