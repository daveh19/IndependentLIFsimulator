#ifndef CL_SYNAPSE_H_
#define CL_SYNAPSE_H_

#include "GeneralIncludes.h"

typedef struct SynapseConstsStruct{
	float gamma_p;
	float gamma_d;
	float theta_p;
	float theta_d;
	unsigned int delay;
	float sigma;
	float tau;
	float tau_ca;
	float c_pre;
	float c_post;
	float dt;
	unsigned int no_syns;
} SynapseConsts;


typedef struct Synapse{
    float * rho;
    float * ca;
	float * gauss;
    //unsigned int ** preT;
	unsigned int * preT;
    unsigned int * postT;
	
	signed int * pre_lif;
	signed int * post_lif;
} cl_Synapse;

typedef struct preT_event_queue{
	int * no_events;
	int ** neuron_id;
} SpikeQueue;


#endif /*CL_SYNAPSE_H_*/
