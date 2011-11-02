#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"
#include "HandleOpenCL.h"


int main (int argc, const char * argv[]) {
	CL cl;
	CL *cl_p = &cl;
	int i;
	
	cl_LIFNeuron lif;
	cl_LIFNeuron *lif_p = &lif;
	(*lif_p).V = malloc(sizeof(float) * DATA_SIZE);
	(*lif_p).I = malloc(sizeof(float) * DATA_SIZE);
	(*lif_p).gauss = calloc(DATA_SIZE, sizeof(float));
	(*lif_p).time_since_spike = calloc(DATA_SIZE, sizeof(unsigned int));
	
	(*lif_p).v_rest = -70.0;
	(*lif_p).v_reset = -68.0;
	(*lif_p).v_threshold = -54.0;
	(*lif_p).r_m = 20.0;
	(*lif_p).c_m = 0.001;
	(*lif_p).sigma = 3.5;
	(*lif_p).refrac_time = 20;
	(*lif_p).dt = 0.001;

	char *KernelSource = readKernelSource("kernel.cl");
	char *k_name = "lif";
	

	printf("initialising data...\n");
//    // Fill our data set with random float values
//    //
//    unsigned int count = DATA_SIZE;
//    for(i = 0; i < count; i++){
//        (*cl_p).data[i] = rand() / (float)RAND_MAX;
//	}
	for ( i = 0; i < DATA_SIZE; i++){
		(*lif_p).V[i] = -70.0;
		(*lif_p).I[i] = 0.3;
		//(*lif_p).gauss = gasdev(&random_seed);
		//(*lif_p).time_since_spike[i] = 0;
	}
	
	
	// OpenCL functions
	if( setupCL(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( makeProgram(cl_p, KernelSource, k_name) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	// OpenCL data IO
	if( createLifIObufs(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( enqueueLifInputBuf(cl_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( setLifKernelArgs(cl_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( getMaxWorkSize(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	// Do the OpenCL processing
	if( enqueueLifKernel(cl_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	waitForKernel(cl_p);
	
	// Read the OpenCL output
	if( enqueueLifOutputBuf(cl_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	//----------------------
	
//    // Validate our results
//    //
//    int correct = 0;
//    for(i = 0; i < count; i++)
//    {
//        if((*cl_p).results[i] == (*cl_p).data[i] * (*cl_p).data[i])
//            correct++;
//		printf("%d %f %f\n", i, (*cl_p).data[i], (*cl_p).results[i]);
//    }
//	
//    // Print a brief summary detailing the results
//    //
//    printf("Computed '%d/%d' correct values!\n", correct, count);	
	
	printf("Results\n");
	for ( i = 0; i < DATA_SIZE; i++){
		printf("V(%d): %f\n", i, (*lif_p).V[i]);
		printf("time_since_spike(%d): %d\n", i, (*lif_p).time_since_spike[i]);
		//(*lif_p).gauss = gasdev(&random_seed);
		//(*lif_p).time_since_spike[i] = 0;
	}
	
	//----------------------
	
	
    shutdownLifKernel(cl_p);
	
    // insert code here...
    printf("Hello, World!\n");
    return 0;
}
