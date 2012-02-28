#ifndef CL_LIFNEURON_H_
#define CL_LIFNEURON_H_

#include "GeneralIncludes.h"

/*
    In order to enable save and resume from a checkpoint export the
    following variables from the command shell:
        LIF_CHECKFILE=filename to save checkpoint information to
        LIF_RESTART=y if this is a resume, unset it otherwise
*/

//int simulation_duration;
//long initial_random_seed;
//long random_seed;
//int resume_offset;
//BOOL checkpointing;
//int siID; // no longer static

//int no_neurons;
//int lif_time_of_last_save;
//int siT; // no longer static
//
//double V_rest;
//double V_reset;
//double V_threshold;
//
//int iRefracTime;
//
//double dCm;
//double dRm;
//double dDt;
//
//double initial_v;
//
//double lifSigma;
//int iPreSpikeDelay;

//typedef int BOOL;

//BOOL a_restart;

//FILE* logfile;
//char* logfilename;
//char logfilearray[FILE_NAME_LENGTH];
//char* neuron_outfilepattern;
//char neuron_outfilearray[FILE_NAME_LENGTH];

//int (*current_fn)(double *, unsigned int);

typedef struct {
    float * V;
	float * I;
	float * gauss;
    unsigned int * time_since_spike;

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

//
//void lif_memory_init(LIFNeuron *);
//int neuron_finalise(int, LIFNeuron *);
//void loadExternalCurrent(LIFNeuron *);
////double calciumFromPreSynapticSpikes(Synapse *);
////double calciumFromPostSynapticSpikes(Synapse *);
////void updateCalciumConcentration(Synapse *);
//// BOOL h(LIFNeuron *, double);
////void updateSynapticEfficacy(Synapse *);
//void updateNeuronMembraneVoltage(LIFNeuron *);

#endif /*CL_LIFNEURON_H_*/
