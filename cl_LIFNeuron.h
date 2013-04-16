#ifndef CL_LIFNEURON_H_
#define CL_LIFNEURON_H_

#include "GeneralIncludes.h"

typedef struct LIFNeuron{
    double * V;
	double * I;
	double * gauss;
    unsigned int * time_since_spike;
	
	unsigned int * no_outgoing_synapses;
	unsigned int * no_outgoing_ee_synapses;
	signed int ** outgoing_synapse_index;
	
	unsigned int * no_incoming_synapses;
	signed int ** incoming_synapse_index;

	double v_rest;
	double v_reset;
	double v_threshold;
	double tau_m;
	//float r_m;
	//float c_m;
	double sigma;
	float refrac_time;  //TODO: why is this a float?
	double dt;
	unsigned int no_lifs;
	
	unsigned int time_step;
	unsigned int random123_seed;
	
	unsigned char * subpopulation_flag; // manipulations will be performed on this population
} cl_LIFNeuron;


#endif /*CL_LIFNEURON_H_*/
