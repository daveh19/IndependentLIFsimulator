/*
 *  DataReporters.c
 *  XclNet
 *
 *  Created by David Higgins on 09/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "DataReporters.h"

void reporters_setup(){
	char outfile[FILE_NAME_LENGTH];
	raster_name = "raster.dat";
	intracellular_name = "intracellular.dat";
	average_activity_name = "network_activity.dat";
	synaptic_activity_name = "single_synapse.dat";
	synaptic_strength_name = "final_synaptic_strength.dat";
	
	// Make sure that directory 'output' exists
	if(mkdir("output",(S_IRUSR | S_IWUSR | S_IXUSR)) == -1){
		if (errno == EEXIST){
			printf("Directory 'output' already exists.\n");
		}
		else{
			perror("Error creating directory 'output'");
		}
	}
	
	// Raster output file
	strcpy(outfile, "output/");
	strcat(outfile, raster_name);
	//printf("DEBUG: %s\n", outfile);
	raster_output = fopen(outfile, "a");
	if(raster_output == NULL){
		perror("Error: failed to open raster output file\n");
	}
	fprintf(raster_output, "\n\n\n\n\n# Raster output (t, lif_no)\n");
	
	// Intracellular recording from a single neuron
	strcpy(outfile, "output/");
	strcat(outfile, intracellular_name);
	//printf("DEBUG: %s\n", outfile);
	intracellular_output = fopen(outfile, "a");
	if(intracellular_output == NULL){
		perror("Error: failed to open intracellular output file\n");
	}
	lif_currents_EE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_EI = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_IE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_II = calloc(MAX_TIME_STEPS, sizeof(float));
	fprintf(intracellular_output, "\n\n\n\n\n# Intracellular recorder (t, V(t), time since spike(t), Iext(t), gauss(t-1), Itot(t)), currents:{EE,IE,EI,II}\n# Neuron ID: %d\n", RECORDER_NEURON_ID);
	
	// Population spiking activity
	strcpy(outfile, "output/");
	strcat(outfile, average_activity_name);
	//printf("DEBUG: %s\n", outfile);
	average_activity_ouput = fopen(outfile, "a");
	if(average_activity_ouput == NULL){
		perror("Error: failed to open average activity output file\n");
	}
	// Setup the bins for recording average population spiking behaviour
	no_spiking_bins = (LIF_DT / BIN_SIZE) * MAX_TIME_STEPS;
	summary_exc_spikes = calloc(no_spiking_bins, sizeof(float));
	summary_inh_spikes = calloc(no_spiking_bins, sizeof(float));
	//TODO: new code for monitoring multiple synapses here
	summary_rho = calloc(no_spiking_bins, sizeof(float));
	summary_M = calloc(no_spiking_bins, sizeof(float));
	summary_S = calloc(no_spiking_bins, sizeof(float));
	summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
	fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin (ms), TotSpikes, ExcSpikes, InhSpikes, InstantaneousExcRate, InstantaneousInhRate, RhoAv, N, RhoStdev)\n# all normalised to their respective population sizes\n");
	
	// Detailed recording from single synapse
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_activity_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_activity_output = fopen(outfile, "a");
	if(synaptic_activity_output == NULL){
		perror("Error: failed to open synaptic activity output file\n");
	}
	fprintf(synaptic_activity_output, "\n\n\n\n\n# Single synapse recorder (t, rho(t), ca(t), preT(t), postT(t))\n# Synapse ID: %d\n", RECORDER_SYNAPSE_ID);
	
	// Final state of all dynamic synapses
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_strength_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_strength_output = fopen(outfile, "a");
	if(synaptic_strength_output == NULL){
		perror("Error: failed to open synaptic strength output file\n");
	}
	fprintf(synaptic_strength_output, "\n\n\n\n\n# Final synaptic strengths (syn_id, pre_syn_lif_id, post_syn_lif_id, rho_initial, rho_final)\n");
	
	#ifdef DEBUG_MODE_NETWORK
		//Debugging stuff
		lif_debug_name = "lif_debug.dat";
		// LIF connectivity and activity output file
		strcpy(outfile, "output/");
		strcat(outfile, lif_debug_name);
		//printf("DEBUG: %s\n", outfile);
		lif_debug_output = fopen(outfile, "a");
		if(lif_debug_output == NULL){
			perror("Error: failed to open lif debug output file\n");
		}
		fprintf(lif_debug_output, "\n\n\n\n\n# LIF debug output (neuron id, mean destination id, no outgoing synapses, sum of gauss, no outgoing EE syns,\n# {no within class, mean dest within class}:{EE,EI,IE,II},\n# no incomming:{EE,EI,IE,II})\n");
		lif_mean_destination = calloc(NO_LIFS, sizeof(float));
		lif_gauss_totals = calloc(NO_LIFS, sizeof(float));
		lif_debug_no_EE = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_EI = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_IE = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_II = calloc(NO_LIFS, sizeof(int));
		lif_mean_dest_EE = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_EI = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_IE = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_II = calloc(NO_LIFS, sizeof(float));
		lif_in_EE = calloc(NO_LIFS, sizeof(int));
		lif_in_EI = calloc(NO_LIFS, sizeof(int));
		lif_in_IE = calloc(NO_LIFS, sizeof(int));
		lif_in_II = calloc(NO_LIFS, sizeof(int));
	#endif /* DEBUG_MODE_NETWORK */
}


void print_raster_spike(int t, int lif_no){
	// A spike has occurred, add its occurrence to raster file
	// t, lif_id
	fprintf(raster_output, "%d %d\n", t, lif_no);
}


void print_network_summary_activity(){
	// bin_id_ms, no_spikes, no_exc_spikes, no_inh_spikes, exc_freq, inh_freq
	printf("Outputting network summary activity\n");
	for(int i = 0; i < no_spiking_bins; i++){
		fprintf(average_activity_ouput, "%d %f %f %f %f %f %f %d %f\n", i, ((summary_inh_spikes[i] + summary_exc_spikes[i]) / NO_LIFS), (summary_exc_spikes[i] / NO_EXC), (summary_inh_spikes[i] / NO_INH), ((summary_exc_spikes[i] / NO_EXC) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), summary_rho[i]/summary_n[i], summary_n[i], sqrt(summary_S[i]/(summary_n[i]-1)) );
	}
}


void print_synapse_activity(int t, cl_Synapse *syn){
	// in event-based model preT and postT here are one timestep in past wrt other variables
	// t, rho, ca, preT, postT
	fprintf(synaptic_activity_output, "%d %f %f %d %d\n", t, (*syn).rho[RECORDER_SYNAPSE_ID], (*syn).ca[RECORDER_SYNAPSE_ID], (*syn).preT[RECORDER_SYNAPSE_ID], (*syn).postT[RECORDER_SYNAPSE_ID]);
}


void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const){
	// syn_id, pre_lif_id, post_lif_id, rho_initial, rho_final
	for(int i = 0; i < (*syn_const).no_syns; i++){
		fprintf(synaptic_strength_output, "%d %d %d %f %f\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho_initial[i], (*syn).rho[i]);
	}
}


void print_lif_debug(cl_LIFNeuron *lif){
	// lif_id, mean_dest_id, no_out_syns, gauss_total, no_out_EE_syns, no_EE, dest_EE, no_EI, dest_EI, no_IE, dest_IE, no_II, dest_II, in_EE, in_EI, in_IE, in_II 
	printf("\nLIF Debug: saving connection statistics to file...\n");
	for(int i = 0; i < (*lif).no_lifs; i++){
		fprintf(lif_debug_output, "%d %f %d %f %d %d %f %d %f %d %f %d %f %d %d %d %d\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i], (*lif).no_outgoing_ee_synapses[i], lif_debug_no_EE[i], lif_mean_dest_EE[i], lif_debug_no_EI[i], lif_mean_dest_EI[i], lif_debug_no_IE[i], lif_mean_dest_IE[i], lif_debug_no_II[i], lif_mean_dest_II[i], lif_in_EE[i], lif_in_EI[i], lif_in_IE[i], lif_in_II[i]);
		//printf("%d %f %d %f\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i]);
	}
}


void reporters_close(){
	fclose(raster_output);
	fclose(intracellular_output);
	fclose(average_activity_ouput);
	fclose(synaptic_activity_output);
	fclose(synaptic_strength_output);
	#ifdef DEBUG_MODE_NETWORK
		fclose(lif_debug_output);
	#endif /* DEBUG_MODE_NETWORK */
}
