#include "GeneralIncludes.h"
#include "HandleOpenCL.h"


int main (int argc, const char * argv[]) {
	CL cl;
	CL *cl_p = &cl;
	int i;

	char *KernelSource = readKernelSource("kernel.cl");
	char *k_name = "square";
	

	printf("loading data...\n");
    // Fill our data set with random float values
    //
    unsigned int count = DATA_SIZE;
    for(i = 0; i < count; i++){
        (*cl_p).data[i] = rand() / (float)RAND_MAX;
	}
	
	
	// OpenCL functions
	if( setupCL(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( makeProgram(cl_p, KernelSource, k_name) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	// OpenCL data IO
	if( createIObufs(cl_p, count) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( enqueueInputBuf(cl_p, count) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( setKernelArgs(cl_p, count) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( getMaxWorkSize(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	// Do the OpenCL processing
	if( enqueueKernel(cl_p, count) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	waitForKernel(cl_p);
	
	// Read the OpenCL output
	if( enqueueOutputBuf(cl_p, count) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	//----------------------
	
    // Validate our results
    //
    int correct = 0;
    for(i = 0; i < count; i++)
    {
        if((*cl_p).results[i] == (*cl_p).data[i] * (*cl_p).data[i])
            correct++;
		printf("%d %f %f\n", i, (*cl_p).data[i], (*cl_p).results[i]);
    }
	
    // Print a brief summary detailing the results
    //
    printf("Computed '%d/%d' correct values!\n", correct, count);	
	
	//----------------------
	
	
    shutdownKernel(cl_p);
	
    // insert code here...
    printf("Hello, World!\n");
    return 0;
}
