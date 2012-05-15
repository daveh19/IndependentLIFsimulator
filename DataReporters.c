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
	printf("DEBUG: %s\n", outfile);
	raster_output = fopen(outfile, "a");
	if(raster_output == NULL){
		perror("Error: failed to open raster output file\n");
	}
	fprintf(raster_output, "\n\n\n\n\n# Raster output (t, lif_no)\n");
	
	// Intracellular recording from a single neuron
	strcpy(outfile, "output/");
	strcat(outfile, intracellular_name);
	printf("DEBUG: %s\n", outfile);
	intracellular_output = fopen(outfile, "a");
	if(intracellular_output == NULL){
		perror("Error: failed to open intracellular output file\n");
	}
	fprintf(intracellular_output, "\n\n\n\n\n# Intracellular recorder (t, V(t), Iext(t), Itot(t))\n# Neuron ID: %d\n", RECORDER_NEURON_ID);
	
	// Population spiking activity
	strcpy(outfile, "output/");
	strcat(outfile, average_activity_name);
	printf("DEBUG: %s\n", outfile);
	average_activity_ouput = fopen(outfile, "a");
	if(average_activity_ouput == NULL){
		perror("Error: failed to open average activity output file\n");
	}
	// Setup the bins for recording average population spiking behaviour
	no_spiking_bins = (LIF_DT / BIN_SIZE) * MAX_TIME_STEPS;
	summary_exc_spikes = calloc(no_spiking_bins, sizeof(float));
	summary_inh_spikes = calloc(no_spiking_bins, sizeof(float));
	fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin (ms), TotSpikes, ExcSpikes, InhSpikes, InstantaneousExcRate, InstantaneousInhRate)\n# all normalised to their respective population sizes\n");
	
	// Detailed recording from single synapse
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_activity_name);
	printf("DEBUG: %s\n", outfile);
	synaptic_activity_output = fopen(outfile, "a");
	if(synaptic_activity_output == NULL){
		perror("Error: failed to open synaptic activity output file\n");
	}
	fprintf(synaptic_activity_output, "\n\n\n\n\n# Single synapse recorder (t, rho(t), ca(t), preT(t), postT(t))\n# Synapse ID: %d\n", RECORDER_SYNAPSE_ID);
	
	// Final state of all dynamic synapses
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_strength_name);
	printf("DEBUG: %s\n", outfile);
	synaptic_strength_output = fopen(outfile, "a");
	if(synaptic_strength_output == NULL){
		perror("Error: failed to open synaptic strength output file\n");
	}
	fprintf(synaptic_strength_output, "\n\n\n\n\n# Final synaptic strengths (syn_id, pre_syn_lif_id, post_syn_lif_id, rho_final)\n");
}


void print_raster_spike(int t, int lif_no){
	// A spike has occurred, add its occurrence to raster file
	fprintf(raster_output, "%d %d\n", t, lif_no);
}


void print_network_summary_activity(){
	printf("Outputting network summary activity\n");
	for(int i = 0; i < no_spiking_bins; i++){
		fprintf(average_activity_ouput, "%d %f %f %f %f %f\n", i, ((summary_inh_spikes[i] + summary_exc_spikes[i]) / NO_LIFS), (summary_exc_spikes[i] / NO_EXC), (summary_inh_spikes[i] / NO_INH), ((summary_exc_spikes[i] / NO_EXC) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)) );
	}
}


void print_synapse_activity(int t, cl_Synapse *syn){
	fprintf(synaptic_activity_output, "%d %f %f %d %d\n", t, (*syn).rho[RECORDER_SYNAPSE_ID], (*syn).ca[RECORDER_SYNAPSE_ID], (*syn).preT[RECORDER_SYNAPSE_ID], (*syn).postT[RECORDER_SYNAPSE_ID]);
}


void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const){
	for(int i = 0; i < (*syn_const).no_syns; i++){
		fprintf(synaptic_strength_output, "%d %d %d %f\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho[i]);
	}
}


void reporters_close(){
	fclose(raster_output);
	fclose(intracellular_output);
	fclose(average_activity_ouput);
	fclose(synaptic_activity_output);
	fclose(synaptic_strength_output);
}
