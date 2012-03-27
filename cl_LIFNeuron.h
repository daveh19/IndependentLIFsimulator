#ifndef CL_LIFNEURON_H_
#define CL_LIFNEURON_H_

#include "GeneralIncludes.h"

typedef struct LIFNeuron{
    float * V;
	float * I;
	float * gauss;
    unsigned int * time_since_spike;
	
	unsigned int * no_outgoing_synapses;
	signed int ** outgoing_synapse_index;
	
	unsigned int * no_incoming_synapses;
	signed int ** incoming_synapse_index;

	float v_rest;
	float v_reset;
	float v_threshold;
	float r_m;
	float c_m;
	float sigma;
	float refrac_time;
	float dt;
	unsigned int no_lifs;
} cl_LIFNeuron;


#endif /*CL_LIFNEURON_H_*/
