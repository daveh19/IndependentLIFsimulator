#ifndef LIFNEURON_H_
#define LIFNEURON_H_

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

int no_neurons;
int lif_time_of_last_save;
int siT; // no longer static

double V_rest;
double V_reset;
double V_threshold;

int iRefracTime;

double dCm;
double dRm;
double dDt;

double initial_v;

double lifSigma;
int iPreSpikeDelay;

//typedef int BOOL;

//BOOL a_restart;

//FILE* logfile;
//char* logfilename;
//char logfilearray[FILE_NAME_LENGTH];
//char* neuron_outfilepattern;
//char neuron_outfilearray[FILE_NAME_LENGTH];

int (*current_fn)(double *, unsigned int);

typedef struct LIFNeuron{
        double * V;
        double * Iext;
        double * Isyn;
        unsigned int * spikeT;
        int ID;
} LIFNeuron;


void lif_memory_init(LIFNeuron *);
int neuron_finalise(int, LIFNeuron *);
void loadExternalCurrent(LIFNeuron *);
//double calciumFromPreSynapticSpikes(Synapse *);
//double calciumFromPostSynapticSpikes(Synapse *);
//void updateCalciumConcentration(Synapse *);
// BOOL h(LIFNeuron *, double);
//void updateSynapticEfficacy(Synapse *);
void updateNeuronMembraneVoltage(LIFNeuron *);

#endif /*LIFNEURON_H_*/
