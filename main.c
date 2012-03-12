#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"
#include "cl_Synapse.h"
#include "HandleOpenCL.h"
#include "NumericalTools.h"


int main (int argc, const char * argv[]) {
	int i, j;
	long random_seed = -13;
	
	// LIF compute kernel
	CL cl_lif;
	CL *cl_lif_p = &cl_lif;
	
	cl_LIFNeuron lif;
	cl_LIFNeuron *lif_p = &lif;
	(*lif_p).V = malloc(sizeof(float) * NO_LIFS);
	(*lif_p).I = malloc(sizeof(float) * NO_LIFS);
	(*lif_p).gauss = calloc(NO_LIFS, sizeof(float));
	(*lif_p).time_since_spike = calloc(NO_LIFS, sizeof(unsigned int));
	
	(*lif_p).v_rest = -70.0;
	(*lif_p).v_reset = -68.0;
	(*lif_p).v_threshold = -54.0;
	(*lif_p).r_m = 20.0;
	(*lif_p).c_m = 0.001;
	(*lif_p).sigma = 0; //3.5;
	(*lif_p).refrac_time = 20;
	(*lif_p).dt = 0.001;
	(*lif_p).no_lifs = NO_LIFS;

	char *KernelSource = readKernelSource("kernel.cl");
	char *k_name_lif = "lif";
	
	// Synapse compute kernel
	CL cl_syn;
	CL *cl_syn_p = & cl_syn;
	
	cl_Synapse syn;
	cl_Synapse *syn_p = &syn;
	(*syn_p).rho = malloc(sizeof(float) * NO_SYNS);
	(*syn_p).ca = malloc(sizeof(float) * NO_SYNS);
	(*syn_p).gauss = calloc(NO_SYNS, sizeof(float));
	//TODO:preT will need to be 2-D to cope with (t-D)
	(*syn_p).preT = calloc(NO_SYNS, sizeof(unsigned int));
	(*syn_p).postT = calloc(NO_SYNS, sizeof(unsigned int));
	
	(*syn_p).gamma_p = 725.085;
	(*syn_p).gamma_d = 331.909;
	(*syn_p).theta_p = 1.3;
	(*syn_p).theta_d = 1.0;
	(*syn_p).delay = 4; // measured in multiples of dt
	(*syn_p).sigma = 3.35;
	(*syn_p).tau = 346.3615;
	(*syn_p).tau_ca = 0.0226936;
	(*syn_p).c_pre = 0.5617539;
	(*syn_p).c_post = 1.23964;
	(*syn_p).dt = 0.001;
	(*syn_p).no_syns = NO_SYNS;
	
	char *k_name_syn = "synapse";
	

	printf("initialising data...\n");
//    // Fill our data set with random float values
//    //
//    unsigned int count = DATA_SIZE;
//    for(i = 0; i < count; i++){
//        (*cl_p).data[i] = rand() / (float)RAND_MAX;
//	}
	for ( i = 0; i < NO_LIFS; i++){
		(*lif_p).V[i] = -66.0;
		(*lif_p).I[i] = 1.0;
		//(*lif_p).gauss[i] = 1.;
		(*lif_p).gauss[i] = gasdev(&random_seed);
		(*lif_p).time_since_spike[i] = (*lif_p).refrac_time;
	}
	for( i = 0; i < NO_SYNS; i++){
		(*syn_p).rho[i] = 1;
		(*syn_p).ca[i] = 5;
		(*syn_p).gauss[i] = gasdev(&random_seed);
	}
	
	
	// OpenCL functions
	if( setupCL(cl_lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( setupCL(cl_syn_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( makeProgram(cl_lif_p, KernelSource, k_name_lif) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( makeProgram(cl_syn_p, KernelSource, k_name_syn) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	// OpenCL data IO
	if( createLifIObufs(cl_lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( createSynIObufs(cl_syn_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	
	if( enqueueLifInputBuf(cl_lif_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( enqueueSynInputBuf(cl_syn_p, syn_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( setLifKernelArgs(cl_lif_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( setSynKernelArgs(cl_syn_p, syn_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( getMaxWorkSize(cl_lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	if( getMaxWorkSize(cl_syn_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}	
	
	// Do the OpenCL processing
	j = 0;
	while(j < MAX_TIME_STEPS){
		// -----Process LIF Kernel-------
		if( enqueueLifKernel(cl_lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}	
		waitForKernel(cl_lif_p);
		// Read the OpenCL output
		if( enqueueLifOutputBuf(cl_lif_p, lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		
		// -----Process Synapse Kernel-----
		if( enqueueSynKernel(cl_syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}	
		waitForKernel(cl_syn_p);
		// Read the OpenCL output
		if( enqueueSynOutputBuf(cl_syn_p, syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}	
	
		// Output results
		printf("V(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[0], j, (*lif_p).time_since_spike[0], (*lif_p).gauss[0]);
		printf("rho(%d): %f, ca(%d): %f, gauss: %f\n", j, (*syn_p).rho[0], j, (*syn_p).ca[0], (*syn_p).gauss[0]);
		
		// Generate new random numbers, for noise processes
		for ( i = 0; i < NO_LIFS; i++){
			(*lif_p).gauss[i] = gasdev(&random_seed);
		}
		for( i = 0; i < NO_SYNS; i++){
			(*syn_p).gauss[i] = gasdev(&random_seed);
		}
		// Setup next LIF Kernel
		if( enqueueLifInputBuf(cl_lif_p, lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		// Setup next Synapse Kernel
		if( enqueueSynInputBuf(cl_syn_p, syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		
		j++;
	}
	//TODO: could do another read here, as bufs have already been enqueued for processing
	
	//----------------------

	/*printf("Results\n");
	for ( i = 0; i < NO_LIFS; i++){
		printf("V(%d): %f\n", i, (*lif_p).V[i]);
		printf("time_since_spike(%d): %d\n", i, (*lif_p).time_since_spike[i]);
		//(*lif_p).gauss = gasdev(&random_seed);
		//(*lif_p).time_since_spike[i] = 0;
	}*/
	/*printf("Results:\n");
	 for ( i = 0; i < NO_LIFS; i++){
	 printf("V(%d): %f, ", i, (*lif_p).V[i]);
	 printf("time_since_spike(%d): %d\n", i, (*lif_p).time_since_spike[i]);
	 }*/
	/*printf("\n --------- \n");
	for ( i = 0; i < NO_SYNS; i++){
		printf("rho(%d): %f, ca(%d): %f\n", i, (*syn_p).rho[i], i, (*syn_p).ca[i]);
	}*/
	
	//----------------------
	
	
    shutdownLifKernel(cl_lif_p);
	shutdownSynKernel(cl_syn_p);
	
    // insert code here...
    printf("Hello, World!\n");
    return 0;
}
