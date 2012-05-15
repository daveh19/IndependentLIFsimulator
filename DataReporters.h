/*
 *  DataReporters.h
 *  XclNet
 *
 *  Created by David Higgins on 09/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralIncludes.h"
#include "cl_Synapse.h"

#define RECORDER_NEURON_ID (0)
#define RECORDER_SYNAPSE_ID (0)


// File pointers for recording data
//Raster output file
char* raster_name;
//Intracellular recorder file
char* intracellular_name;
//Averaged population activity recorder file
char* average_activity_name;// = "av_activity.dat";
//Individual synaptic activity file
char* synaptic_activity_name;// = "synapse.dat";
//
char* synaptic_strength_name;// = "synapse_summary.dat";

FILE *raster_output;
FILE *intracellular_output;
FILE *average_activity_ouput;
FILE *synaptic_activity_output;
FILE *synaptic_strength_output;

float *summary_exc_spikes;
float *summary_inh_spikes;
unsigned int no_spiking_bins;

void reporters_setup();
void reporters_close();

void print_raster_spike(int t, int lif_no);
void print_network_summary_activity();
void print_synapse_activity(int t, cl_Synapse *syn);
void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const);