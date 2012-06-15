
//	(*syn_p).gamma_p = 725.085;
//	(*syn_p).gamma_d = 331.909;
//	(*syn_p).theta_p = 1.3;
//	(*syn_p).theta_d = 1.0;
//	//(*syn_p).delay = 4; // measured in multiples of dt
//	(*syn_p).sigma = 0; //3.35; //TODO: switch noise back on
//	(*syn_p).tau = 346.3615;
//	(*syn_p).tau_ca = 0.0226936;
//	(*syn_p).c_pre = 0.5617539;
//	(*syn_p).c_post = 1.23964;
//	(*syn_p).dt = 0.001;
//	//(*syn_p).no_syns = NO_SYNS;

/*
 // Generate links within Synapses
 // CONSIDER: this code could be much faster by using a temporarily much larger memory space and doing
 // processing withing generation of LIF links loops, then copying or resizing to real memory requirements
 printf("...links from synapses to pre and post synaptic neurons\n");
 (*syn_p).pre_lif = calloc(total_synapses, sizeof(signed int));
 (*syn_p).post_lif = calloc(total_synapses, sizeof(signed int));
 network_seed = NETWORK_SEED;
 int k = 0;
 for(int i = 0; i < NO_LIFS; i++){
 for(int j = 0; j < NO_LIFS; j++){
 if(i != j){
 if ((ran2(&network_seed)) < p){
 (*syn_p).pre_lif[k] = i;
 (*syn_p).post_lif[k] = j;
 k++;
 //printf("(%d)->(%d)\n", i, j);
 }
 }
 }
 }
 printf("k final (total synapses): %d\n", k);
 */	


//printf("Number of outgoing connections from lif(%d): %d\n", i, (*lif_p).no_outgoing_synapses[i]);
/*(*lif_p).outgoing_synapse_index[i] = calloc((*lif_p).no_outgoing_synapses[i], sizeof(signed int));
 for(int j = 0; j < k; j++){
 (*lif_p).outgoing_synapse_index[i][j] = out_links[j];
 //printf("Index(%d)(%d): %d \n", i, j, links[j]);
 }*/







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
#ifdef DEBUG_MODE
	printf("V(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[RECORDER_NEURON_ID], j, (*lif_p).time_since_spike[RECORDER_NEURON_ID], (*lif_p).gauss[RECORDER_NEURON_ID]);
	printf("rho(%d): %f, ca(%d): %f, preT(%d): %d, postT(%d): %d, gauss: %f\n", j, (*syn_p).rho[RECORDER_NEURON_ID], j, (*syn_p).ca[RECORDER_NEURON_ID], j, (*syn_p).preT[RECORDER_NEURON_ID], j, (*syn_p).postT[RECORDER_NEURON_ID], (*syn_p).gauss[RECORDER_NEURON_ID]);
#endif /* DEBUG_MODE */
	
	// ---- Prepare next run ----
	
	/*
	 //TODO: Event-based 1 (Delayed forward propagation)
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
	 //TODO: Event-based 2 (Reset delayed event queue)
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
			//TODO: Event-based 3 (Backpropagation to synapse)
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
				//TODO: Event-based 4 (Update in advance of current transfer)
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
			//TODO: Event-based 5 (Add spike to delayed processing event queue)
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
		// Alternative to event queue system, assumes only 1 spike can occur in delay period
		// Only one of these two systems can be used at a time, currently using event queue system
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
	fprintf(intracellular_output, "%f %f %f %f %f\n", (*lif_p).I[RECORDER_NEURON_ID], lif_currents_EE[j], lif_currents_IE[j], lif_currents_EI[j], lif_currents_II[j]);
	
	// Print state of a single synapse
	print_synapse_activity(j, syn_p);
	
	
#ifdef DEBUG_MODE
	printf("after transfer V(%d): %f, I(%d): %f, time_since_spike(%d): %d, gauss: %f\n", j, (*lif_p).V[RECORDER_NEURON_ID], j, (*lif_p).I[RECORDER_NEURON_ID], j, (*lif_p).time_since_spike[RECORDER_NEURON_ID], (*lif_p).gauss[RECORDER_NEURON_ID]);
	printf("after transfer rho(%d): %f, ca(%d): %f, preT(%d): %d, postT(%d): %d, gauss: %f\n", j, (*syn_p).rho[RECORDER_NEURON_ID], j, (*syn_p).ca[RECORDER_NEURON_ID], j, (*syn_p).preT[RECORDER_NEURON_ID], j, (*syn_p).postT[RECORDER_NEURON_ID], (*syn_p).gauss[RECORDER_NEURON_ID]);
#endif /* DEBUG_MODE */
	
	// Setup next LIF Kernel
	// this is the part I was able to comment out and sim still worked! (most of the time!)
	if( enqueueLifInputBuf(cl_lif_p, lif_p, rnd_lif_p) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	/*
	 // Setup next Synapse Kernel
	 if( enqueueSynInputBuf(cl_syn_p, syn_p, syn_const_p, rnd_syn_p) == EXIT_FAILURE){
	 return EXIT_FAILURE;
	 }
	 */
	
	//TODO: Event-based 6 (Update event queue offset variable)
	//offset = (++offset) % (*syn_const_p).delay;
	j++;
}