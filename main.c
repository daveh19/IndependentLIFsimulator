#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"
#include "cl_Synapse.h"
#include "HandleOpenCL.h"
#include "NumericalTools.h"

#include "DataReporters.h"

//#define PI (atan(1.)*4.)

void freeMemory(cl_LIFNeuron *lif_p, cl_Synapse *syn_p, FixedSynapse *fixed_syn_p, SpikeQueue *spike_queue_p);


unsigned int generateNetwork(cl_LIFNeuron *lif_p, cl_Synapse *syn_p, FixedSynapse *fixed_syn_p){
	//clock_t start, finish;
	//double totaltime;
	float p = CONNECTIVITY_PROBABILITY;
	long network_seed = NETWORK_SEED;
	unsigned int total_ee_synapses = 0;
	unsigned int total_fixed_synapses = 0;
	
	unsigned int mean_ee_synapses_per_neuron = (int)((NO_EXC)*(CONNECTIVITY_PROBABILITY));
	unsigned int mean_total_synapses_per_neuron = (int)((NO_LIFS)*(CONNECTIVITY_PROBABILITY));
	//CONSIDER: widening the margin of error on the following two estimates
	unsigned int estimated_total_ee_synapses = (NO_EXC * (mean_ee_synapses_per_neuron+100)) + (int)(NO_EXC / 10) + 100; // mean + wide margin + constant (for small nets)
	unsigned int estimated_total_synapses = (NO_LIFS * (mean_total_synapses_per_neuron+100)) + (int)(NO_LIFS / 10) + 100; // mean + wide margin + constant (for small nets)
	unsigned int estimated_ee_synapses_per_neuron = (mean_ee_synapses_per_neuron) + (int)(mean_ee_synapses_per_neuron/2) + 1000; // mean + wide margin + constant (for small nets)
	unsigned int estimated_total_synapses_per_neuron = (mean_total_synapses_per_neuron) + (int)(mean_total_synapses_per_neuron/2) + 2000; // mean + wide margin + constant (for small nets)
	
	float delta_spike_modifier = ((*lif_p).r_m * (*lif_p).c_m) / (*lif_p).dt;
	//printf("DEBUG: delta_spike_modifier %f\n", delta_spike_modifier);
	
	(*syn_p).pre_lif = calloc(estimated_total_ee_synapses, sizeof(signed int));
	(*syn_p).post_lif = calloc(estimated_total_ee_synapses, sizeof(signed int));
	
	(*fixed_syn_p).Jx = calloc(estimated_total_synapses, sizeof(float));
	(*fixed_syn_p).post_lif = calloc(estimated_total_synapses, sizeof(signed int));
	
	printf("mean_ee_syn_per_neuron: %d, est_ee_syn_per_neuron: %d, est_total_ee_synapses: %d, est_total_syn_per_neuron: %d\n", mean_ee_synapses_per_neuron, estimated_ee_synapses_per_neuron, estimated_total_ee_synapses, estimated_total_synapses_per_neuron);
	
	printf("Generating network structure...\n");
	//start = clock();
	
	// Assign basic memory requirements for keeping track of pre and post neuronal synapses
	(*lif_p).no_outgoing_synapses = calloc(NO_LIFS, sizeof(unsigned int));
	(*lif_p).no_outgoing_ee_synapses = calloc(NO_LIFS, sizeof(unsigned int));
	(*lif_p).outgoing_synapse_index = malloc(sizeof(signed int *) * NO_LIFS);
	
	(*lif_p).no_incoming_synapses = calloc(NO_LIFS, sizeof(unsigned int));
	(*lif_p).incoming_synapse_index = malloc(sizeof(signed int *) * NO_LIFS);
	
	// Assign (hopefully) overly large memory for recording ids of pre and post neuronal synapses
	for(int i = 0; i < NO_LIFS; i++){
		(*lif_p).incoming_synapse_index[i] = malloc(sizeof(signed int) * estimated_total_synapses_per_neuron);
		(*lif_p).outgoing_synapse_index[i] = malloc(sizeof(signed int) * estimated_total_synapses_per_neuron);
	}
	
	// Generate Synapses randomly, telling each synapse its pre and post synaptic neuron ids and each lif its pre or post neuronal synapse ids
	for(int i = 0; i < NO_EXC; i++){
		// By generating EE synapses first we can assume that all synapses with array address < total_ee_synapses are EE synapses
		for(int j = 0; j < NO_EXC; j++){
			// EXC -> EXC synapses
			if(i != j){ // Disallow autapses
				if ((ran2(&network_seed)) < p){ 
					// A new synapse is formed
					
					//printf("synapse(%d) ", total_synapses);
					// Assign indices of pre and post synaptic neurons to the new synapse
					(*syn_p).pre_lif[total_ee_synapses] = i;
					(*syn_p).post_lif[total_ee_synapses] = j;
					//printf("pre_lif: %d, post_lif: %d, ", (*syn_p).pre_lif[total_synapses], (*syn_p).post_lif[total_synapses]);
					
					// Update pre-synaptic neuron relationship with synapse
					(*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]] = (int)total_ee_synapses; //synaptic id
					(*lif_p).no_outgoing_synapses[i]++;// could be added to array lookup in previous line
					(*lif_p).no_outgoing_ee_synapses[i]++;
					//printf("out_id: %d, no_out: %d ", (*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]-1], (*lif_p).no_outgoing_synapses[i]);
					
					// Update post-synaptic neuron relationship with synapse					
					(*lif_p).incoming_synapse_index[j][(*lif_p).no_incoming_synapses[j]] = (int)total_ee_synapses; //syn id
					(*lif_p).no_incoming_synapses[j]++;
					//printf("in_id: %d, no_in: %d \n", (*lif_p).incoming_synapse_index[j][(*lif_p).no_incoming_synapses[j]-1], (*lif_p).no_incoming_synapses[j]);
					
					total_ee_synapses++;
					//Debug code
					/*lif_mean_destination[i] += j;
					lif_mean_dest_EE[i] += j;
					lif_debug_no_EE[i]++;
					lif_in_EE[j]++;*/
				}
			}
		}
	}
	for(int i = 0; i < NO_EXC; i++){
		for(int j = NO_EXC; j < NO_LIFS; j++){
			// EXC -> INH synapses
			if(i != j){ //TODO: test not strictly necessary here
				if((ran2(&network_seed)) < p){
					// A new synapse
					// by not having a pre_lif variable we save on checks when backpropagating spikes
					(*fixed_syn_p).post_lif[total_fixed_synapses] = j;
					(*fixed_syn_p).Jx[total_fixed_synapses] = delta_spike_modifier * J_IE;
				
					(*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]] = (int)(total_fixed_synapses + total_ee_synapses); //synaptic id
					(*lif_p).no_outgoing_synapses[i]++;
				
					total_fixed_synapses++;
					//Debug code
					/*lif_mean_destination[i] += j;
					lif_mean_dest_IE[i] += j;
					lif_debug_no_IE[i]++;
					lif_in_IE[j]++;*/
				}
			}
		}
	}
	for(int i = NO_EXC; i < NO_LIFS; i++){
		for(int j = 0; j < NO_EXC; j++){
			// INH -> EXC synapses
			if(i != j){ //TODO: test not strictly necessary here
				if((ran2(&network_seed)) < p){
					// A new synapse
					(*fixed_syn_p).post_lif[total_fixed_synapses] = j;
					(*fixed_syn_p).Jx[total_fixed_synapses] = delta_spike_modifier * J_EI;
					
					(*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]] = (int)(total_fixed_synapses + total_ee_synapses); //synaptic id
					(*lif_p).no_outgoing_synapses[i]++;
					
					total_fixed_synapses++;
					//Debug code
					/*lif_mean_destination[i] += j;
					lif_mean_dest_EI[i] += j;
					lif_debug_no_EI[i]++;
					lif_in_EI[j]++;*/
				}
			}
		}
		for(int j = NO_EXC; j < NO_LIFS; j++){
			// INH -> INH synapses
			if(i != j){
				if((ran2(&network_seed)) < p){
					// A new synapse
					(*fixed_syn_p).post_lif[total_fixed_synapses] = j;
					(*fixed_syn_p).Jx[total_fixed_synapses] = delta_spike_modifier * J_II;
					
					(*lif_p).outgoing_synapse_index[i][(*lif_p).no_outgoing_synapses[i]] = (int)(total_fixed_synapses + total_ee_synapses); //synaptic id
					(*lif_p).no_outgoing_synapses[i]++;
					
					total_fixed_synapses++;	
					//Debug code
					/*lif_mean_destination[i] += j;
					lif_mean_dest_II[i] += j;
					lif_debug_no_II[i]++;
					lif_in_II[j]++;*/
				}
			}
		}
	}

	//finish = clock();
	//totaltime = (double)(finish - start)/CLOCKS_PER_SEC;
	//printf("Time taken: %f, ", totaltime);
	printf("Total EE synapses: %d, total fixed synapses: %d\n", total_ee_synapses, total_fixed_synapses);
	
	
	// Shrink memory reserved to required sizes
	(*syn_p).pre_lif = realloc((*syn_p).pre_lif, sizeof(signed int) * total_ee_synapses);
	(*syn_p).post_lif = realloc((*syn_p).post_lif, sizeof(signed int) * total_ee_synapses);
	(*fixed_syn_p).Jx = realloc((*fixed_syn_p).Jx, sizeof(float) * total_fixed_synapses);
	(*fixed_syn_p).post_lif = realloc((*fixed_syn_p).post_lif, sizeof(signed int) * total_fixed_synapses);
	//FILE *fp = fopen("output/test_network.dat", "a");
	//fprintf(fp, "# Test of net (no incomming EE, no outgoing, no outgoing EE, (no out - no out EE))\n");
	for(int i = 0; i < NO_LIFS; i++){
		(*lif_p).incoming_synapse_index[i] = realloc((*lif_p).incoming_synapse_index[i], sizeof(signed int) * (*lif_p).no_incoming_synapses[i]);
		(*lif_p).outgoing_synapse_index[i] = realloc((*lif_p).outgoing_synapse_index[i], sizeof(signed int) * (*lif_p).no_outgoing_synapses[i]);
		//Debug
		//fprintf(fp, "%d %d %d %d %d\n", i, (*lif_p).no_incoming_synapses[i], (*lif_p).no_outgoing_synapses[i], (*lif_p).no_outgoing_ee_synapses[i], ((*lif_p).no_outgoing_synapses[i] - (*lif_p).no_outgoing_ee_synapses[i]));
		/*lif_mean_destination[i] /= (*lif_p).no_outgoing_synapses[i];
		lif_mean_dest_EE[i] /= lif_debug_no_EE[i];
		lif_mean_dest_IE[i] /= lif_debug_no_IE[i];
		lif_mean_dest_EI[i] /= lif_debug_no_EI[i];
		lif_mean_dest_II[i] /= lif_debug_no_II[i];*/
	}
	//fclose(fp);
	
	(*fixed_syn_p).total_fixed_synapses = total_fixed_synapses;
	
	return total_ee_synapses;
}


int main (int argc, const char * argv[]) {
	int i, j, k;
	int offset;
	
	clock_t start_t,finish_t;
	double totaltime;

	//Setup output files
	reporters_setup();
	
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
	
	(*lif_p).v_rest = LIF_V_REST;
	(*lif_p).v_reset = LIF_V_RESET;
	(*lif_p).v_threshold = LIF_V_THRESHOLD;
	(*lif_p).r_m = LIF_RM;
	(*lif_p).c_m = LIF_CM;
	(*lif_p).sigma = LIF_SIGMA; //5; 
	(*lif_p).refrac_time = LIF_REFRAC_TIME; //20;
	(*lif_p).dt = LIF_DT;
	(*lif_p).no_lifs = NO_LIFS;
	(*cl_lif_p).job_size = (*lif_p).no_lifs;
	
	// Setup external and synaptic voltages/currents
	float external_voltage = J_EXT;
	// Syanptic currents must be modified by (tau_m/dt) as they are delta current spikes
	float delta_spike_modifier = ((*lif_p).r_m * (*lif_p).c_m) / (*lif_p).dt;
	float transfer_voltage = J_EE;
	transfer_voltage *= delta_spike_modifier;
	//printf("DEBUG: delta_spike_modifier %f, transfer_voltage %f\n", delta_spike_modifier, transfer_voltage);
	 
	char *k_name_lif = "lif";
	
	
	// Synapse compute kernel
	CL cl_syn;
	CL *cl_syn_p = &cl_syn;
	cl_Synapse syn;
	cl_Synapse *syn_p = &syn;
	SynapseConsts syn_const;
	SynapseConsts *syn_const_p = &syn_const;
	
	// Fixed strength excitatory and inhibitory synapses
	FixedSynapse fixed_syn;
	FixedSynapse *fixed_syn_p = &fixed_syn;
	
	// Network generation
	start_t = clock();
	(*syn_const_p).no_syns = generateNetwork(lif_p, syn_p, fixed_syn_p);
	finish_t = clock();
	totaltime = (double)(finish_t - start_t)/CLOCKS_PER_SEC;
	//(*syn_const_p).no_syns = NO_SYNS;
	(*cl_syn_p).job_size = (*syn_const_p).no_syns;
	printf("Network generated, jobsize: %d, generation time: %f\n", (*cl_syn_p).job_size, totaltime);	
	//DEBUGGING code: print network description using different mapping approaches
	// Traverse pre-synaptic neurons
	/*for(int i = 0; i < NO_LIFS; i++){
		for(int j = 0; j < (*lif_p).no_outgoing_ee_synapses[i]; j++){
			printf("LIF(%d) -> LIF(%d), via Synapse(%d), no outgoing %d, no EE outgoing: %d\n", i, (*syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][j]], (*lif_p).outgoing_synapse_index[i][j], (*lif_p).no_outgoing_synapses[i], (*lif_p).no_outgoing_ee_synapses[i]);
		}
		for(int j = (*lif_p).no_outgoing_ee_synapses[i]; j < (*lif_p).no_outgoing_synapses[i]; j++){
			printf("LIF(%d) -> LIF(%d), via fixed Synapse(%d), no outgoing %d, no EE outgoing: %d\n", i, (*fixed_syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][j]-(*syn_const_p).no_syns], ((*lif_p).outgoing_synapse_index[i][j] - (*syn_const_p).no_syns), (*lif_p).no_outgoing_synapses[i], (*lif_p).no_outgoing_ee_synapses[i]);
		}
	} */ /*
	printf("---------------\n");
	// Traverse post-synaptic neurons
	for(int i = 0; i < NO_LIFS; i++){
		for(int j = 0; j < (*lif_p).no_incoming_synapses[i]; j++){
			printf("LIF(%d) <- LIF(%d), via Synapse(%d), no incoming %d\n", i, (*syn_p).pre_lif[(*lif_p).incoming_synapse_index[i][j]], (*lif_p).incoming_synapse_index[i][j], (*lif_p).no_incoming_synapses[i]);
			
		}
	} */ /*
	printf("---------------\n");
	// Traverse synapses
	for(int i = 0; i < (*syn_const_p).no_syns; i++){
		printf("Synapse(%d) links LIF(%d) -> LIF(%d)\n", i, (*syn_p).pre_lif[i], (*syn_p).post_lif[i]);
	}
	for(int i = 0; i < (*fixed_syn_p).total_fixed_synapses; i++){
		printf("Fixed Synapse(%d) links to LIF(%d)\n", i, (*fixed_syn_p).post_lif[i]);
	}*/
	//End DEBUGGING code
	//End network generation
	
	(*syn_p).rho = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).ca = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).gauss = calloc((*syn_const_p).no_syns, sizeof(float));
	(*syn_const_p).delay = SYN_CALCIUM_DELAY; // measured in multiples of dt

	(*syn_p).preT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	(*syn_p).postT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	
	// Event queue for delayed propagation of pre-synaptic spikes to synaptic calcium buffer
	SpikeQueue spike_queue;
	SpikeQueue *spike_queue_p = &spike_queue;
	(*spike_queue_p).neuron_id = malloc((*syn_const_p).delay * sizeof(int *));
	(*spike_queue_p).no_events = calloc((*syn_const_p).delay, sizeof(int));
	for(i = 0; i < (*syn_const_p).delay; i++){
		(*spike_queue_p).neuron_id[i] = malloc((*lif_p).no_lifs * sizeof(int));
	}
	

	(*syn_const_p).gamma_p = SYN_GAMMA_P;
	(*syn_const_p).gamma_d = SYN_GAMMA_D;
	(*syn_const_p).theta_p = SYN_THETA_P;
	(*syn_const_p).theta_d = SYN_THETA_D;
	(*syn_const_p).sigma = SYN_SIGMA;
	(*syn_const_p).tau = SYN_TAU;
	(*syn_const_p).tau_ca = SYN_TAU_CA;
	(*syn_const_p).c_pre = SYN_C_PRE;
	(*syn_const_p).c_post = SYN_C_POST;
	(*syn_const_p).dt = SYN_DT;
	
	
	char *k_name_syn = "synapse";
	
	// Random number generator streams
	//for LIF
	cl_MarsagliaStruct rnd_lif;
	cl_MarsagliaStruct *rnd_lif_p = &rnd_lif;
	(*rnd_lif_p).d_z = malloc((*lif_p).no_lifs * sizeof(unsigned int));
	(*rnd_lif_p).d_w = malloc((*lif_p).no_lifs * sizeof(unsigned int));
	(*rnd_lif_p).d_jsr = malloc((*lif_p).no_lifs * sizeof(unsigned int));
	(*rnd_lif_p).d_jcong = malloc((*lif_p).no_lifs * sizeof(unsigned int));
	//for Synapse
	cl_MarsagliaStruct rnd_syn;
	cl_MarsagliaStruct *rnd_syn_p = &rnd_syn;
	(*rnd_syn_p).d_z = malloc((*syn_const_p).no_syns * sizeof(unsigned int));
	(*rnd_syn_p).d_w = malloc((*syn_const_p).no_syns * sizeof(unsigned int));
	(*rnd_syn_p).d_jsr = malloc((*syn_const_p).no_syns * sizeof(unsigned int));
	(*rnd_syn_p).d_jcong = malloc((*syn_const_p).no_syns * sizeof(unsigned int));
	
	
	printf("initialising data...\n");
    // Prepopulate data set, including with random values
    //
	for ( i = 0; i < (*lif_p).no_lifs; i++){
		//CONSIDER: initialising V and time_since_spike to random values (within reasonable ranges)
		(*lif_p).V[i] = LIF_V_INITIAL;
		(*lif_p).I[i] = external_voltage;
		(*lif_p).time_since_spike[i] = (*lif_p).refrac_time;
		(*rnd_lif_p).d_z[i] = 362436069 + i + 1 + PARALLEL_SEED;
		(*rnd_lif_p).d_w[i] = 521288629 + i + 1 + PARALLEL_SEED;
		(*rnd_lif_p).d_jsr[i] = 123456789 + i + 1 + PARALLEL_SEED;
		(*rnd_lif_p).d_jcong[i] = 380116160 + i + 1 + PARALLEL_SEED;
	}
	for( i = 0; i < (*syn_const_p).no_syns; i++){
		//TODO: store rho_init, when it becomes a non constant initial value
		(*syn_p).rho[i] = SYN_RHO_INITIAL;
		(*syn_p).ca[i] = SYN_CA_INITIAL;
		(*rnd_syn_p).d_z[i] = 362436069 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_w[i] = 521288629 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_jsr[i] = 123456789 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_jcong[i] = 380116160 - i + PARALLEL_SEED;
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
	
	
	if( enqueueLifInputBuf(cl_lif_p, lif_p, rnd_lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( enqueueSynInputBuf(cl_syn_p, syn_p, syn_const_p, rnd_syn_p) == EXIT_FAILURE){
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
	offset = 0;
	//clock_t start,finish;
	//double totaltime;
	printf("Go\n");
	start_t = clock();
	while(j < MAX_TIME_STEPS){
		
		// -----Process LIF Kernel-------
		if( enqueueLifKernel(cl_lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		//TODO: reenable waitForKernel() ?
		//waitForKernel(cl_lif_p);
		// Read the OpenCL output
		if( enqueueLifOutputBuf(cl_lif_p, lif_p, rnd_lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		
		/*
		// -----Process Synapse Kernel-----
		if( enqueueSynKernel(cl_syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}	
		waitForKernel(cl_syn_p);
		// Read the OpenCL output
		if( enqueueSynOutputBuf(cl_syn_p, syn_p, syn_const_p, rnd_syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		 */
		 
	
		// Output results
		//TODO: commented for debug
		//printf("V(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[RECORDER_NEURON_ID], j, (*lif_p).time_since_spike[RECORDER_NEURON_ID], (*lif_p).gauss[RECORDER_NEURON_ID]);
		//printf("rho(%d): %f, ca(%d): %f, preT(%d): %d, postT(%d): %d, gauss: %f\n", j, (*syn_p).rho[RECORDER_NEURON_ID], j, (*syn_p).ca[RECORDER_NEURON_ID], j, (*syn_p).preT[RECORDER_NEURON_ID], j, (*syn_p).postT[RECORDER_NEURON_ID], (*syn_p).gauss[RECORDER_NEURON_ID]);
		//(*syn_p).preT[0][0]
		
		// ---- Prepare next run ----
		
		/*
		// Transfer delayed pre-synaptic spikes to EE synapses
		for( i = 0; i < (*spike_queue_p).no_events[offset]; i++){
			//printf("number of events: %d\n", (*spike_queue_p).no_events[offset]);
			//printf("No outgoing EE synapses for neuron(%d): %d\n", (*spike_queue_p).neuron_id[offset][i], (*lif_p).no_outgoing_ee_synapses[(*spike_queue_p).neuron_id[offset][i]]);
			
			// Process each neuron which spiked (delay timesteps ago)
			for ( k = 0; k < (*lif_p).no_outgoing_ee_synapses[ (*spike_queue_p).neuron_id[offset][i] ]; k++){
				//if(i==0){
				//	printf("transferring delayed pre-synaptic spike to synapse(%d)\n", (*lif_p).outgoing_synapse_index[ (*spike_queue_p).neuron_id[offset][i] ][k]); 
				//}
				(*syn_p).preT[ (*lif_p).outgoing_synapse_index[ (*spike_queue_p).neuron_id[offset][i] ][k] ] = 1;
			}
		}
		// Reset delayed event queue
		(*spike_queue_p).no_events[offset] = 0;
		 */

		// Apply external voltage (this cannot be done in same loop as spike detection/propagation
		for( i = 0; i < (*lif_p).no_lifs; i++){
			// Fixed external current
			(*lif_p).I[i] = external_voltage;
			//Debug code
			//lif_gauss_totals[i] += (*lif_p).gauss[i];
		}
		
		// Print to intracellular recorder file
		// print: time, voltage, input current
		fprintf(intracellular_output, "%d %f %d %f %f ", j, (*lif_p).V[RECORDER_NEURON_ID], (*lif_p).time_since_spike[RECORDER_NEURON_ID], (*lif_p).I[RECORDER_NEURON_ID], (*lif_p).gauss[RECORDER_NEURON_ID]);
		
		//int local_count = 0;
		// Update LIFs: spike detection/propagation to post-synaptic lifs as well as pre- and post-lif neurons
		for ( i = 0; i < (*lif_p).no_lifs; i++){
			if((*lif_p).time_since_spike[i] == 0){
				//CONSIDER: using local variables to point to I[], post_lif[], Jx[], etc. it cuts down on dereferencing!
				//TODO: strongly consider implementing parallel spike transfer system
				/*if(i==0){
					printf("%d spiked\n", i);
				}*/
				// Post-synaptic spike should backpropagate to its synapses with no delay
				/*for ( k = 0; k < (*lif_p).no_incoming_synapses[i]; k++){
					// as non EE based lifs are not added to incoming_synapses lists this is safe
					(*syn_p).postT[(*lif_p).incoming_synapse_index[i][k]] = 1;
					//if(i==0){
					//	printf("backprop to pre-lif synapse(%d)\n", (*lif_p).incoming_synapse_index[i][k]);
					//}
				}*/
				// Transfer voltage change to post-synaptic neurons
				for ( k = 0; k < (*lif_p).no_outgoing_ee_synapses[i]; k++){
					// across plastic synapses
					//Debug code
					if ((*syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][k]] == RECORDER_NEURON_ID){
						//local_count++;
						lif_currents_EE[j] += transfer_voltage * (*syn_p).rho[(*lif_p).outgoing_synapse_index[i][k]]; 
					}
					(*lif_p).I[(*syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][k]]] += transfer_voltage * (*syn_p).rho[(*lif_p).outgoing_synapse_index[i][k]]; 
					/*if(i==0){
						printf("current transfer, I: %f, to post-synaptic neuron(%d)\n", (transfer_voltage * (*syn_p).rho[(*lif_p).outgoing_synapse_index[i][k]]), (*syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][k]]);
					}*/				
				}
				for ( k = (*lif_p).no_outgoing_ee_synapses[i]; k < (*lif_p).no_outgoing_synapses[i]; k++){
					// across fixed synapses
					//Debug code
					if((*fixed_syn_p).post_lif[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ] == RECORDER_NEURON_ID){
						//local_count--;
						if((i < NO_EXC) && (RECORDER_NEURON_ID >= NO_EXC)){ //E->I
							lif_currents_IE[j] += (*fixed_syn_p).Jx[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ];
						}
						else if((i >= NO_EXC) && (RECORDER_NEURON_ID >= NO_EXC)){ //I->I
							lif_currents_II[j] += (*fixed_syn_p).Jx[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ];
						}
						else if((i >= NO_EXC) && (RECORDER_NEURON_ID < NO_EXC)){ //I->E
							lif_currents_EI[j] += (*fixed_syn_p).Jx[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ];
						}
					}
					// Alternative dereferencing for fixed synapses
					(*lif_p).I[(*fixed_syn_p).post_lif[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ] ] += (*fixed_syn_p).Jx[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ];
					//if(i==0){
					//	printf("current transfer, I: %f, via fixed syn(%d) to post-synaptic neuron(%d)\n", ((*fixed_syn_p).Jx[ ((*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns) ]), (*lif_p).outgoing_synapse_index[i][k]-(*syn_const_p).no_syns, (*fixed_syn_p).post_lif[(*lif_p).outgoing_synapse_index[i][k] - (*syn_const_p).no_syns]);
					//}					
				}
				// Add to pre-spike event queue
				//CONSIDER: don't add non EE events to event queue (relative efficiencies depend on NO_INH<<NO_EXC and nu_i>nu_e)
				/*(*spike_queue_p).neuron_id[offset][(*spike_queue_p).no_events[offset]] = i;
				(*spike_queue_p).no_events[offset]++;
				 */
				
				//Print to raster file
				print_raster_spike(j, i);
				
				// Add to average spiking activity bins
				if(i < NO_EXC){
					summary_exc_spikes[(int)( ( (*lif_p).dt / BIN_SIZE ) * j + EPSILLON)]++;
				}
				else{
					summary_inh_spikes[(int)( ( (*lif_p).dt / BIN_SIZE ) * j + EPSILLON)]++;
				}
			} // end of handling spike
			// Pre-synaptic spike propagates across synapse after delay
			//Alternative to event queue system, assumes only 1 spike can occur in delay period
			/*else if((*lif_p).time_since_spike[i] == (*syn_const_p).delay){
				for ( k = 0; k < (*lif_p).no_outgoing_synapses[i]; k++){
					(*syn_p).preT[(*lif_p).outgoing_synapse_index[i][k]] = 1;
					if(i==0){
						printf("Delayed spike transferred to post-lif synapse(%d)\n", (*lif_p).outgoing_synapse_index[i][k]);
					}
				}
			}
			*/
		} // end of loop over neurons
		//printf("count: %d\n", local_count);
		
		// Print total I to intracellular recorder file
		//if(i == RECORDER_NEURON_ID){
		fprintf(intracellular_output, "%f %f %f %f %f\n", (*lif_p).I[RECORDER_NEURON_ID], lif_currents_EE[j], lif_currents_IE[j], lif_currents_EI[j], lif_currents_II[j]);
		//}
		
		// Print state of a single synapse
		//print_synapse_activity(j, syn_p);
		
		
		//TODO: commented for debug
		//printf("after transfer V(%d): %f, I(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[RECORDER_NEURON_ID], j, (*lif_p).I[RECORDER_NEURON_ID], j, (*lif_p).time_since_spike[RECORDER_NEURON_ID], (*lif_p).gauss[RECORDER_NEURON_ID]);
		//printf("after transfer rho(%d): %f, ca(%d): %f, preT(%d): %d, postT(%d): %d, gauss: %f\n", j, (*syn_p).rho[RECORDER_NEURON_ID], j, (*syn_p).ca[RECORDER_NEURON_ID], j, (*syn_p).preT[RECORDER_NEURON_ID], j, (*syn_p).postT[RECORDER_NEURON_ID], (*syn_p).gauss[RECORDER_NEURON_ID]);
		
		
		// Setup next LIF Kernel
		// TODO: this is the part I am able to comment out and sim still works! (most of the time!)
		if( enqueueLifInputBuf(cl_lif_p, lif_p, rnd_lif_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		/*
		// Setup next Synapse Kernel
		if( enqueueSynInputBuf(cl_syn_p, syn_p, syn_const_p, rnd_syn_p) == EXIT_FAILURE){
			return EXIT_FAILURE;
		}
		 */
		
		//offset = (++offset) % (*syn_const_p).delay;
		j++;
	}
	finish_t = clock();
	printf("Stop\n");
	totaltime = (double)(finish_t - start_t)/CLOCKS_PER_SEC;
	printf("Main loop run time: %lf secs, start time: %lf, finish time: %lf\n", totaltime, (double)start_t, (double)finish_t);
	
	//TODO: could do another process and read here, as bufs have already been enqueued for processing
	
	printf("Simulation finished, printing summary of network activity...\n");
	// Print summary of excitatory and inhibitory activity
	print_network_summary_activity();
	
	printf("done.\nAnd final state of synapses...");
	// Print final state of synapse strengths
	//print_synapses_final_state(syn_p, syn_const_p);
	//printf("done.\n");
	
	//Debug
	//print_lif_debug(lif_p);
	
	printf("done\n");
	// Close output files
	reporters_close();
	
    shutdownLifKernel(cl_lif_p);
	shutdownSynKernel(cl_syn_p);
	
	freeMemory(lif_p, syn_p, fixed_syn_p, spike_queue_p);
	
    printf("Hello, World!\n");
    return 0;
}


void freeMemory(cl_LIFNeuron *lif_p, cl_Synapse *syn_p, FixedSynapse *fixed_syn_p, SpikeQueue *spike_queue_p){
	// LIF variables
	free((*lif_p).V);
	free((*lif_p).I);
	free((*lif_p).gauss);
	free((*lif_p).time_since_spike);
	free((*lif_p).no_outgoing_synapses);
	free((*lif_p).no_outgoing_ee_synapses);
	free((*lif_p).no_incoming_synapses);
	for(int i = 0; i < (*lif_p).no_lifs; i++){
		free((*lif_p).outgoing_synapse_index[i]);
		free((*lif_p).incoming_synapse_index[i]);
	}
	free((*lif_p).outgoing_synapse_index);
	free((*lif_p).incoming_synapse_index);
	
	// Synapse variables
	free((*syn_p).rho);
	free((*syn_p).ca);
	free((*syn_p).gauss);
	free((*syn_p).preT);
	free((*syn_p).postT);
	free((*syn_p).pre_lif);
	free((*syn_p).post_lif);
	
	free((*fixed_syn_p).Jx);
	free((*fixed_syn_p).post_lif);
	
	free((*spike_queue_p).no_events);
	for(int i = 0; i < SYN_CALCIUM_DELAY; i++){
		free((*spike_queue_p).neuron_id[i]);
	}
	free((*spike_queue_p).neuron_id);
	
	// Reporter variables
	free(summary_exc_spikes);
	free(summary_inh_spikes);
	
	// Debugging variables
	free(lif_gauss_totals);
	free(lif_mean_destination);
	free(lif_debug_no_EE);
	free(lif_debug_no_IE);
	free(lif_debug_no_EI);
	free(lif_debug_no_II);
	free(lif_mean_dest_EE);
	free(lif_mean_dest_IE);
	free(lif_mean_dest_EI);
	free(lif_mean_dest_II);
	free(lif_in_EE);
	free(lif_in_IE);
	free(lif_in_EI);
	free(lif_in_II);
	free(lif_currents_EE);
	free(lif_currents_IE);
	free(lif_currents_EI);
	free(lif_currents_II);
}
