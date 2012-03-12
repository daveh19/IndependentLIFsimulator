#ifndef CL_SYNAPSE_H_
#define CL_SYNAPSE_H_

#include "GeneralIncludes.h"

/*
    In order to enable save and resume from a checkpoint export the
    following variables from the command shell:
        SYNAPSE_CHECKFILE=filename to save checkpoint information to
        SYNAPSE_RESTART=y if this is a resume, unset it otherwise
*/

/*
int simulation_duration;
int no_synapses;

double initial_c;
double initial_rho;

long initial_random_seed;
long random_seed;

int iTau;
int iTauC;

double dRhoFixed;

double dCpre;
double dCpost;
double dThetaD;
double dThetaP;
double dGammaD; //CONSIDER: does this really need to be double?
double dGammaP; //CONSIDER: does this really need to be double?

double dSigma;
int iPreSpikeDelay;
double poisson_param;

int siT; // no longer static
int siID; // no longer static
int time_of_last_save;
int resume_offset;

typedef int BOOL;

BOOL checkpointing;
//BOOL a_restart;

FILE* logfile;
char* logfilename;
char* outfilepattern;
char logfilearray[FILE_NAME_LENGTH];
char outfilearray[FILE_NAME_LENGTH];

int (*train_fn)(unsigned int *, unsigned int *, unsigned int);
*/

typedef struct Synapse{
    float * rho;
    float * ca;
	float * gauss;
    unsigned int ** preT;
    unsigned int * postT;
	
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
} cl_Synapse;


/*
void synapse_memory_init(Synapse *);
int finalise(int, Synapse *);
void loadInitialSpikeTimes(Synapse *);
double calciumFromPreSynapticSpikes(Synapse *);
double calciumFromPostSynapticSpikes(Synapse *);
void updateCalciumConcentration(Synapse *);
BOOL h(Synapse *, double);
void updateSynapticEfficacy(Synapse *);
 */

#endif /*CL_SYNAPSE_H_*/
