#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"
#include "cl_Synapse.h"
#include "HandleOpenCL.h"
#include "NumericalTools.h"

//#define PI (atan(1.)*4.)

void freeMemory(cl_LIFNeuron *lif_p, cl_Synapse *syn_p);

unsigned int generateNetwork(cl_LIFNeuron *lif_p, cl_Synapse *syn_p){
	float p = CONNECTIVITY_PROBABILITY;
	long network_seed = NETWORK_SEED;
	unsigned int total_synapses = 0;
	
	unsigned int mean_synapses_per_neuron = (int)((NO_LIFS)*(CONNECTIVITY_PROBABILITY));
	//CONSIDER: widening the margin of error on the following two estimates
	unsigned int estimated_total_synapses = (NO_LIFS * mean_synapses_per_neuron) + (int)(NO_LIFS / 10) + 100; // mean + wide margin + constant (for small nets)
	unsigned int estimated_synapses_per_neuron = (mean_synapses_per_neuron) + (int)(mean_synapses_per_neuron/2) + 100; // mean + wide margin + constant (for small nets)
	
	(*syn_p).pre_lif = calloc(estimated_total_synapses, sizeof(signed int));
	(*syn_p).post_lif = calloc(estimated_total_synapses, sizeof(signed int));
	
	printf("mean_syn_per_neuron: %d, est_syn_per_neuron: %d, est_total_synapses: %d\n", mean_synapses_per_neuron, estimated_synapses_per_neuron, estimated_total_synapses);
	
	printf("Generating network structure...\n");
	
	// Assign basic memory requirements for keeping track of pre and post neuronal synapses
	(*lif_p).no_outgoing_synapses = calloc(NO_LIFS, sizeof(unsigned int));
	(*lif_p).outgoing_synapse_index = malloc(sizeof(signed int *) * NO_LIFS);
	
	(*lif_p).no_incoming_synapses = calloc(NO_LIFS, sizeof(unsigned int));
	(*lif_p).incoming_synapse_index = malloc(sizeof(signed int *) * NO_LIFS);
	
	// Assign (hopefully) overly large memory for recording ids of pre and post neuronal synapses
	for(int i = 0; i < NO_LIFS; i++){
		(*lif_p).incoming_synapse_index[i] = malloc(sizeof(signed int) * estimated_synapses_per_neuron);
		(*lif_p).outgoing_synapse_index[i] = malloc(sizeof(signed int) * estimated_synapses_per_neuron);
	}
	
	// Generate Synapses randomly, telling each synapse its pre and post synaptic neuron ids and each lif its pre or post neuronal synapse ids
	for(int i = 0; i < NO_LIFS; i++){
		for(int j = 0; j < NO_LIFS; j++){
			if(i != j){
				if ((ran2(&network_seed)) < p){ 
					// A new synapse is formed
					
					//printf("synapse(%d) ", total_synapses);
					// Assign indices of pre and post synaptic neurons to the new synapse
					(*syn_p).pre_lif[total_synapses] = i;
					(*syn_p).post_lif[total_synapses] = j;
					//printf("pre_lif: %d, post_lif: %d, ", (*syn_p).pre_lif[total_synapses], (*syn_p).post_lif[total_synapses]);
					
					// Update pre-synaptic neuron relationship with synapse
					(*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]] = (int)total_synapses; //synaptic id
					(*lif_p).no_outgoing_synapses[i]++;// could be added to array lookup in previous line
					//printf("out_id: %d, no_out: %d ", (*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]-1], (*lif_p).no_outgoing_synapses[i]);
					
					// Update post-synaptic neuron relationship with synapse					
					(*lif_p).incoming_synapse_index[j][(*lif_p).no_incoming_synapses[j]] = (int)total_synapses; //syn id
					(*lif_p).no_incoming_synapses[j]++;
					//printf("in_id: %d, no_in: %d \n", (*lif_p).incoming_synapse_index[j][(*lif_p).no_incoming_synapses[j]-1], (*lif_p).no_incoming_synapses[j]);
					
					total_synapses++;
				}
			}
		}

	}

	printf("Total synapses: %d\n", total_synapses);
	
	
	// Shrink memory reserved to required sizes
	(*syn_p).pre_lif = realloc((*syn_p).pre_lif, sizeof(signed int) * total_synapses);
	(*syn_p).post_lif = realloc((*syn_p).post_lif, sizeof(signed int) * total_synapses);
	for(int i = 0; i < NO_LIFS; i++){
		(*lif_p).incoming_synapse_index[i] = realloc((*lif_p).incoming_synapse_index[i], sizeof(signed int) * (*lif_p).no_incoming_synapses[i]);
		(*lif_p).outgoing_synapse_index[i] = realloc((*lif_p).outgoing_synapse_index[i], sizeof(signed int) * (*lif_p).no_incoming_synapses[i]);
	}
	
	return total_synapses;
}

int main (int argc, const char * argv[]) {
	int i, j;
	long random_seed = -13;

	char *KernelSource = readKernelSource("kernel.cl");
	
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
	(*lif_p).sigma = 0; //3.5; //TODO: switch noise back on
	(*lif_p).refrac_time = 20;
	(*lif_p).dt = 0.001;
	(*lif_p).no_lifs = NO_LIFS;
	(*cl_lif_p).job_size = (*lif_p).no_lifs;
	 
	char *k_name_lif = "lif";
	
	
	// Synapse compute kernel
	CL cl_syn;
	CL *cl_syn_p = &cl_syn;
	cl_Synapse syn;
	cl_Synapse *syn_p = &syn;
	SynapseConsts syn_const;
	SynapseConsts *syn_const_p = &syn_const;
	
	// Network generation
	(*syn_const_p).no_syns = generateNetwork(lif_p, syn_p);
	//(*syn_const_p).no_syns = NO_SYNS;
	(*cl_syn_p).job_size = (*syn_const_p).no_syns;
	printf("Network generated, jobsize: %d\n", (*cl_syn_p).job_size);
	/*for(int i = 0; i < NO_LIFS; i++){
		for(int j = 0; j < (*lif_p).no_outgoing_synapses[i]; j++){
			printf("LIF(%d) -> LIF(%d), via Synapse(%d), expected no outgoing %d\n", i, (*syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][j]], (*lif_p).outgoing_synapse_index[i][j], (*lif_p).no_outgoing_synapses[i]);
		}
	}
	printf("---------------\n");
	for(int i = 0; i < NO_LIFS; i++){
		for(int j = 0; j < (*lif_p).no_incoming_synapses[i]; j++){
			printf("LIF(%d) <- LIF(%d), via Synapse(%d), expected no incoming %d\n", i, (*syn_p).pre_lif[(*lif_p).incoming_synapse_index[i][j]], (*lif_p).incoming_synapse_index[i][j], (*lif_p).no_incoming_synapses[i]);
			
		}
	}
	printf("---------------\n");
	for(int i = 0; i < (*syn_const_p).no_syns; i++){
		printf("Synapse(%d) links LIF(%d) -> LIF(%d)\n", i, (*syn_p).pre_lif[i], (*syn_p).post_lif[i]);
	}*/
	//End network generation
	
	(*syn_p).rho = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).ca = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).gauss = calloc((*syn_const_p).no_syns, sizeof(float));
	(*syn_const_p).delay = CALCIUM_DELAY; // measured in multiples of dt

	/*(*syn_p).preT = malloc(sizeof(unsigned int *) * (*syn_const_p).delay);
	for(i = 1; i < (*syn_const_p).delay; i++){
		(*syn_p).preT[i] = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	}*/
	(*syn_p).preT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	(*syn_p).postT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	

	(*syn_const_p).gamma_p = 725.085;
	(*syn_const_p).gamma_d = 331.909;
	(*syn_const_p).theta_p = 1.3;
	(*syn_const_p).theta_d = 1.0;
	(*syn_const_p).sigma = 0; //3.35; //TODO: switch noise back on
	(*syn_const_p).tau = 346.3615;
	(*syn_const_p).tau_ca = 0.0226936;
	(*syn_const_p).c_pre = 0.5617539;
	(*syn_const_p).c_post = 1.23964;
	(*syn_const_p).dt = 0.001;
	
	
	char *k_name_syn = "synapse";
	

	printf("initialising data...\n");
    // Prepopulate data set, including with random values
    //
	for ( i = 0; i < NO_LIFS; i++){
		(*lif_p).V[i] = -66.0;
		(*lif_p).I[i] = 1.0;
		//(*lif_p).gauss[i] = 1.;
		(*lif_p).gauss[i] = gasdev(&random_seed);
		(*lif_p).time_since_spike[i] = (*lif_p).refrac_time;
	}
	for( i = 0; i < (*syn_const_p).no_syns; i++){
		(*syn_p).rho[i] = 1;
		(*syn_p).ca[i] = 5;
		(*syn_p).gauss[i] = gasdev(&random_seed);
		/*for( j = 0; j < (*syn_const_p).delay; j++){
			(*syn_p).preT[j][i] = j;//0;//j;
		}*/
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
	if( enqueueSynInputBuf(cl_syn_p, syn_p, syn_const_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( setLifKernelArgs(cl_lif_p, lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( setSynKernelArgs(cl_syn_p, syn_p, syn_const_p) == EXIT_FAILURE){
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
		if( enqueueSynOutputBuf(cl_syn_p, syn_p, syn_const_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		 
	
		// Output results
		printf("V(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[0], j, (*lif_p).time_since_spike[0], (*lif_p).gauss[0]);
		printf("rho(%d): %f, ca(%d): %f, preT(%d): %d, postT(%d): %d, gauss: %f\n", j, (*syn_p).rho[0], j, (*syn_p).ca[0], j, (*syn_p).preT[0], j, (*syn_p).postT[0], (*syn_p).gauss[0]);
		//(*syn_p).preT[0][0]
		
		//TODO: when Network is in place, change preT and postT to pointers to time_since_spike in relevant lif's
		/*printf("BEFORE:-------------\n");
		for( i = 0; i < (*syn_p).delay; i++){
			for(int k = 0; k < 10; k++){
				printf("preT[%d][%d]: %d, ", i, k, (*syn_p).preT[i][k]);
			}
			printf("\n");
		}*/
		
		//Rearrange pointers in preT[] to take account of delays
		/*unsigned int * local_delay_rearrange = (*syn_p).preT[0];
		for( i = 0; i < ((*syn_const_p).delay - 1); i++){
			(*syn_p).preT[i] = (*syn_p).preT[i+1];
		}
		(*syn_p).preT[(*syn_const_p).delay - 1] = local_delay_rearrange;
		*/
		
		/*printf("AFTER:-------------\n");
		for( i = 0; i < (*syn_p).delay; i++){
			for(int k = 0; k < 10; k++){
				printf("preT[%d][%d]: %d, ", i, k, (*syn_p).preT[i][k]);
			}
			printf("\n");
		}*/
		
		//TODO: hook up corresponding time_since_spike vars to synapse preT and postT
		
		// Generate new random numbers, for noise processes
		for ( i = 0; i < (*lif_p).no_lifs; i++){
			(*lif_p).gauss[i] = gasdev(&random_seed);
		}
		for( i = 0; i < (*syn_const_p).no_syns; i++){
			(*syn_p).gauss[i] = gasdev(&random_seed);
		}
		
		// Setup next LIF Kernel
		if( enqueueLifInputBuf(cl_lif_p, lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		// Setup next Synapse Kernel
		if( enqueueSynInputBuf(cl_syn_p, syn_p, syn_const_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		
		j++;
	}
	//TODO: could do another process and read here, as bufs have already been enqueued for processing
	
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
	
	freeMemory(lif_p, syn_p);
	
    printf("Hello, World!\n");
    return 0;
}

void freeMemory(cl_LIFNeuron *lif_p, cl_Synapse *syn_p){
	free((*lif_p).V);
	free((*lif_p).I);
	free((*lif_p).gauss);
	free((*lif_p).time_since_spike);
	free((*lif_p).no_outgoing_synapses);
	free((*lif_p).no_incoming_synapses);
	for(int i = 0; i < (*lif_p).no_lifs; i++){
		free((*lif_p).outgoing_synapse_index[i]);
		free((*lif_p).incoming_synapse_index[i]);
	}
	free((*lif_p).outgoing_synapse_index);
	free((*lif_p).incoming_synapse_index);
	
	free((*syn_p).rho);
	free((*syn_p).ca);
	free((*syn_p).gauss);
	free((*syn_p).preT);
	free((*syn_p).postT);
	free((*syn_p).pre_lif);
	free((*syn_p).post_lif);
}
