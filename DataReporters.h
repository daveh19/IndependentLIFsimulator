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

//Debug
#include "cl_LIFNeuron.h"


// File pointers for recording data
// Raster output file
char* raster_name; //"raster.dat"
// Intracellular recorder file
char* intracellular_name; //"intracellular.dat"
// Averaged population activity file
char* average_activity_name;// = "network_activity.dat";
// Individual synaptic activity file
char* synaptic_activity_name;// = "single_synapse.dat";
// Final synaptic strengths of all dynamic synapses file
char* synaptic_strength_name;// = "final_synaptic_strength.dat";

// Calculation of synchange for each frequency
char* synchange_name;// = "synchange.dat";

//TODO: rename file pointer variable names for DataReporters
FILE *raster_output;
FILE *intracellular_output;
FILE *average_activity_ouput; //network_activity_output
FILE *synaptic_activity_output; //single_synapse_output
FILE *synaptic_strength_output; //final_synaptic_strength_output

FILE *synchange_output; //synchange_output

// Summary variables for monitoring network firing rate
//CONSIDER: since we use a timestepping approach these variables could be condensed
// to single value variables and printed out during the simulation
float *summary_exc_spikes;
//float *summary_inh_spikes;
unsigned int no_spiking_bins;

// Summary variables for monitoring multiple recorder synapses
/*float *summary_rho;
float *summary_M;
float *summary_S;
unsigned int *summary_n;*/
// Summary variables for main population synapse recorders
float *pop_summary_rho;
float *pop_summary_M;
float *pop_summary_S;
unsigned int *pop_summary_n;

// Variables for manipulating subset of neurons
//float *lif_injection_spikes;
//int no_injection_lifs;

//Debugging variables
float *lif_gauss_totals;
//float *lif_mean_destination;
//char* lif_debug_name;
//FILE *lif_debug_output;
/*int *lif_debug_no_EE;
int *lif_debug_no_IE;
int *lif_debug_no_EI;
int *lif_debug_no_II;
float *lif_mean_dest_EE;
float *lif_mean_dest_IE;
float *lif_mean_dest_EI;
float *lif_mean_dest_II;
int *lif_in_EE;
int *lif_in_IE;
int *lif_in_EI;
int *lif_in_II;
float *lif_currents_EE;
float *lif_currents_IE;
float *lif_currents_EI;
float *lif_currents_II;*/

void reporters_setup();
void reporters_close();
void reporters_flush();
void alloc_reporter_variables();
void free_reporter_variables();

void print_synchange(cl_Synapse *syn, SynapseConsts *syn_const, double fup, double cmich, double nT);

void print_raster_spike(int t, int lif_no, float isi);
void print_network_summary_activity();
void print_synapse_activity(int t, cl_Synapse *syn);
void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const);
//void print_lif_debug(cl_LIFNeuron *lif);