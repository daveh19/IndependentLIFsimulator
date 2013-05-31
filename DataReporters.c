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
	
	synchange_name = "synchange.dat";
	
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
	fprintf(raster_output, "\n\n\n\n\n# Raster output (t, lif_no, isi)\n");
	
	// Intracellular recording from a single neuron
	strcpy(outfile, "output/");
	strcat(outfile, intracellular_name);
	//printf("DEBUG: %s\n", outfile);
	intracellular_output = fopen(outfile, "a");
	if(intracellular_output == NULL){
		perror("Error: failed to open intracellular output file\n");
	}
	/*lif_currents_EE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_EI = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_IE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_II = calloc(MAX_TIME_STEPS, sizeof(float));*/
	fprintf(intracellular_output, "\n\n\n\n\n# Intracellular recorder (t, V(t), time since spike(t), Iext(t), gauss(t-1), Itot(t))\n# Neuron ID: %d\n", RECORDER_NEURON_ID);
	
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
	//summary_exc_spikes = calloc(no_spiking_bins, sizeof(float)); //Moved to alloc_reporter_variables()
	//summary_inh_spikes = calloc(no_spiking_bins, sizeof(float));
	// new code for monitoring multiple synapses here
	/*summary_rho = calloc(no_spiking_bins, sizeof(float));
	summary_M = calloc(no_spiking_bins, sizeof(float));
	summary_S = calloc(no_spiking_bins, sizeof(float));
	summary_n = calloc(no_spiking_bins, sizeof(unsigned int));*/
	// new code for monitoring main population synapses here
	// moved following to allocate_reporter_variables()
	/*pop_summary_rho = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_rho == NULL){
		printf("Failed to allocate memory for pop_summary_rho\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_M = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_M == NULL){
		printf("Failed to allocate memory for pop_summary_M\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_S = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_S == NULL){
		printf("Failed to allocate memory for pop_summary_S\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
	if(pop_summary_n == NULL){
		printf("Failed to allocate memory for pop_summary_n\n");
		exit(EXIT_FAILURE);
	}*/
	// Code for monitoring selectively manipulated neurons
	//lif_injection_spikes = calloc(no_spiking_bins, sizeof(float));
	fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin, ExcSpikes/N_EXC, InstantaneousExcRate, RhoAv, RhoStdev, NoExcSpikes)\n# all normalised to their respective population sizes\n");
	// Detailed recording from single synapse
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_activity_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_activity_output = fopen(outfile, "a");
	if(synaptic_activity_output == NULL){
		perror("Error: failed to open synaptic activity output file\n");
	}
	fprintf(synaptic_activity_output, "\n\n\n\n\n# Single synapse recorder (t, rho(t), ca(t), preT(t), postT(t))\n# Synapse ID: %d\n", RECORDER_SYNAPSE_ID);
	
	printf("DEBUG: flushing single_synapse.dat\n");
	fflush(synaptic_activity_output);
	
	// Final state of all dynamic synapses
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_strength_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_strength_output = fopen(outfile, "a");
	if(synaptic_strength_output == NULL){
		perror("Error: failed to open synaptic strength output file\n");
	}
	fprintf(synaptic_strength_output, "\n\n\n\n\n# Final synaptic strengths (syn_id, pre_syn_lif_id, post_syn_lif_id, rho_initial, rho_final, alpha_d, alpha_p)\n");
	
	/*#ifdef DEBUG_MODE_NETWORK
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
		fprintf(lif_debug_output, "\n\n\n\n\n# LIF debug output (neuron id, mean destination id, no outgoing synapses, sum of gauss, no outgoing EE syns,\n# {no within class, mean dest within class}:{EE,EI,IE,II},\n# no incomming:{EE,EI,IE,II}), in_sub_pop\n");
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
	#endif*/ /* DEBUG_MODE_NETWORK */
	
	// Synaptic change calculation output
	strcpy(outfile, "output/");
	strcat(outfile, synchange_name);
	//printf("DEBUG: %s\n", outfile);
	synchange_output = fopen(outfile, "a");
	if (synchange_output == NULL){
		perror("Error: failed to open synchange output file\n");
	}
	fprintf(synchange_output, "\n\n\n\n\n# Synchange output (rate, alpha_d, alpha_p, synchange, rhobar)\n");	
}


void print_synchange(cl_Synapse *syn, SynapseConsts *syn_const, double fup, double cmich, double nT){
	double mean_alpha_d = 0;
	double mean_alpha_p = 0;
	double mean_rate = 0;
	
	// Calculate mean values for alpha_d and alpha_p
	for(int i = 0; i < (*syn_const).no_syns; i++){
		//TODO: temporary modification to discard first 3 secs
		mean_alpha_d += (*syn).alpha_d[i] / ((*syn).time_of_last_update[i]*(*syn_const).dt - 3);
		mean_alpha_p += (*syn).alpha_p[i] / ((*syn).time_of_last_update[i]*(*syn_const).dt - 3);
	}
	mean_alpha_d /= (double)(*syn_const).no_syns;
	mean_alpha_p /= (double)(*syn_const).no_syns;
	
	
	// Calculate mean firing rate from 3rd time bin to end
	for(int i = 3; i < no_spiking_bins; i++){
		mean_rate += (summary_exc_spikes[i] / NO_EXC) * (1. / BIN_SIZE);
	}
	mean_rate /= (double)(no_spiking_bins-3);
	
	
	// Begin synchange calculation code
	double rhobar, sigmap, taueff, UP, DOWN, synchange;
	
	double rhostar = 0.5;
	
	rhobar = ((*syn_const).gamma_p * mean_alpha_p) / (((*syn_const).gamma_p * mean_alpha_p) + ((*syn_const).gamma_d * mean_alpha_d));
	
	sigmap = (*syn_const).sigma * sqrt( (mean_alpha_p + mean_alpha_d) / (((*syn_const).gamma_p * mean_alpha_p) + ((*syn_const).gamma_d * mean_alpha_d)) );
	taueff = (*syn_const).tau / (((*syn_const).gamma_p * mean_alpha_p) + ((*syn_const).gamma_d * mean_alpha_d));
	UP = 0.5*(1-erf((rhostar-rhobar+rhobar*exp(-nT/(mean_rate*taueff)))/(sigmap*sqrt(1-exp(-(2*nT)/(mean_rate*taueff))))));
	DOWN = 0.5*(1+erf((rhostar-rhobar+(rhobar-1)*exp(-nT/(mean_rate*taueff)))/(sigmap*sqrt(1.-exp(-(2*nT)/(mean_rate*taueff))))));
	synchange = ( (fup * (1-DOWN) + (1-fup) * UP) * cmich + fup * DOWN + (1-fup) * (1-UP) ) / (fup * cmich + 1 - fup);
	
	fprintf(synchange_output, "%f %f %f %f %f\n", mean_rate, mean_alpha_d, mean_alpha_p, synchange, rhobar);
}


//TODO: add access to (*lif_p) and print out ISI
void print_raster_spike(int t, int lif_no, float isi){
	// A spike has occurred, add its occurrence to raster file
	// print inter-spike-interval too
	// t, lif_id, isi
	fprintf(raster_output, "%d %d %f\n", t, lif_no, isi);
}


void print_network_summary_activity(){
	// bin_id_ms, no_spikes, no_exc_spikes, no_inh_spikes, exc_freq, inh_freq
	printf("Outputting network summary activity\n");
	for(int i = 0; i < no_spiking_bins; i++){
		//fprintf(average_activity_ouput, "%d %f %f %f %f %f %f %f %f %f %f %d\n", i, ((summary_inh_spikes[i] + summary_exc_spikes[i]) / NO_LIFS), (summary_exc_spikes[i] / NO_EXC), (summary_inh_spikes[i] / NO_INH), ((summary_exc_spikes[i] / NO_EXC) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), summary_rho[i]/summary_n[i], sqrt(summary_S[i]/(summary_n[i]-1)),  pop_summary_rho[i]/pop_summary_n[i], sqrt(pop_summary_S[i]/(pop_summary_n[i]-1)), ((lif_injection_spikes[i] / no_injection_lifs) * (1.0 / BIN_SIZE)), (int) summary_exc_spikes[i] );
		fprintf(average_activity_ouput, "%d %f %f %f %f %d\n", i, (summary_exc_spikes[i] / NO_EXC), ((summary_exc_spikes[i] / NO_EXC) * (1.0 / BIN_SIZE)),  pop_summary_rho[i]/pop_summary_n[i], sqrt(pop_summary_S[i]/(pop_summary_n[i]-1)), (int) summary_exc_spikes[i] );
	}
	fprintf(average_activity_ouput, "\n\n\n\n\n");
}


void print_synapse_activity(int t, cl_Synapse *syn){
	// in event-based model preT and postT here are one timestep in past wrt other variables
	// t, rho, ca, preT, postT
	fprintf(synaptic_activity_output, "%d %f %f %d %d\n", t, (*syn).rho[RECORDER_SYNAPSE_ID], (*syn).ca[RECORDER_SYNAPSE_ID], (*syn).preT[RECORDER_SYNAPSE_ID], (*syn).postT[RECORDER_SYNAPSE_ID]);

	//TODO: consider adding flush() to output here
	//fflush(synaptic_activity_output);
}


/*void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const){
	// syn_id, pre_lif_id, post_lif_id, rho_initial, rho_final
	for(int i = 0; i < (*syn_const).no_syns; i++){
		fprintf(synaptic_strength_output, "%d %d %d %f %f\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho_initial[i], (*syn).rho[i]);
	}
}*/

void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const){
	// syn_id, pre_lif_id, post_lif_id, rho_initial, rho_final
	for(int i = 0; i < (*syn_const).no_syns; i++){
		fprintf(synaptic_strength_output, "%d %d %d %f %f %.8f %.8f\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho_initial[i], (*syn).rho[i], (*syn).alpha_d[i]/((*syn).time_of_last_update[i]*(*syn_const).dt), (*syn).alpha_p[i]/((*syn).time_of_last_update[i]*(*syn_const).dt));
	}
	fprintf(synaptic_strength_output, "\n\n\n\n\n");
}

void reporters_flush(){
	fprintf(synaptic_activity_output, "\n\n\n\n\n");
	fprintf(intracellular_output, "\n\n\n\n\n");
	fprintf(raster_output, "\n\n\n\n\n");
	
	fflush(stdout);
	
	fflush(raster_output);
	fflush(intracellular_output);
	fflush(average_activity_ouput);
	fflush(synaptic_activity_output);
	fflush(synaptic_strength_output);
	fflush(synchange_output);
}

/*void print_lif_debug(cl_LIFNeuron *lif){
	// lif_id, mean_dest_id, no_out_syns, gauss_total, no_out_EE_syns, no_EE, dest_EE, no_EI, dest_EI, no_IE, dest_IE, no_II, dest_II, in_EE, in_EI, in_IE, in_II 
	printf("\nLIF Debug: saving connection statistics to file...\n");
	for(int i = 0; i < (*lif).no_lifs; i++){
		fprintf(lif_debug_output, "%d %f %d %f %d %d %f %d %f %d %f %d %f %d %d %d %d %d\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i], (*lif).no_outgoing_ee_synapses[i], lif_debug_no_EE[i], lif_mean_dest_EE[i], lif_debug_no_EI[i], lif_mean_dest_EI[i], lif_debug_no_IE[i], lif_mean_dest_IE[i], lif_debug_no_II[i], lif_mean_dest_II[i], lif_in_EE[i], lif_in_EI[i], lif_in_IE[i], lif_in_II[i], (*lif).subpopulation_flag[i]);
		//printf("%d %f %d %f\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i]);
	}
}*/

void alloc_reporter_variables(){
	//no_spiking_bins = (LIF_DT / BIN_SIZE) * MAX_TIME_STEPS;
	summary_exc_spikes = calloc(no_spiking_bins, sizeof(float));
	if(summary_exc_spikes == NULL){
		printf("Failed to allocate memory for summary_exc_spikes\n");
		exit(EXIT_FAILURE);
	}
	
	pop_summary_rho = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_rho == NULL){
		printf("Failed to allocate memory for pop_summary_rho\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_M = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_M == NULL){
		printf("Failed to allocate memory for pop_summary_M\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_S = calloc(no_spiking_bins, sizeof(float));
	if(pop_summary_S == NULL){
		printf("Failed to allocate memory for pop_summary_S\n");
		exit(EXIT_FAILURE);
	}
	pop_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
	if(pop_summary_n == NULL){
		printf("Failed to allocate memory for pop_summary_n\n");
		exit(EXIT_FAILURE);
	}
	
	/*summary_rho = calloc(no_spiking_bins, sizeof(float));
	 summary_M = calloc(no_spiking_bins, sizeof(float));
	 summary_S = calloc(no_spiking_bins, sizeof(float));
	 summary_n = calloc(no_spiking_bins, sizeof(unsigned int));*/
}

void free_reporter_variables(){
	free(summary_exc_spikes);
	
	free(pop_summary_rho);
	free(pop_summary_M);
	free(pop_summary_S);
	free(pop_summary_n);
	
	/*free(summary_inh_spikes);
	 free(lif_currents_EE);
	 free(lif_currents_IE);
	 free(lif_currents_EI);
	 free(lif_currents_II);*/
	
	/*#ifdef DEBUG_MODE_NETWORK
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
	#endif*/ /* DEBUG_MODE_NETWORK */	
}

void reporters_close(){
	fclose(raster_output);
	fclose(intracellular_output);
	fclose(average_activity_ouput);
	fclose(synaptic_activity_output);
	fclose(synaptic_strength_output);
	fclose(synchange_output);
	/*#ifdef DEBUG_MODE_NETWORK
		fclose(lif_debug_output);
	#endif*/ /* DEBUG_MODE_NETWORK */
}
