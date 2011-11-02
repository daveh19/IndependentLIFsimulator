#ifndef DATATOOLS_H_
#define DATATOOLS_H_

#import "Synapse.h"

int saveSynapseProgressToFile(char* filename, void *obj, int end_time);
int createOutputFileHeader(char* filename, void *obj, int duration, double dCpre, double dCpost, double dThetaD, double dThetaP, double dGammaD, double dGammaP, double dSigma, int iPreSpikeDelay, int iTau, int iTauC, double dRhoFixed, double poisson_param, long initial_random_seed);
int checkpoint_save(Synapse *syn);

int checkpoint_load(FILE *checkpoint_fp, Synapse *syn);

Synapse* checkpoint_init(int argc, char *argv[], Synapse *syn);

void loadSimulationParameters(int argc, char *argv[]);
int printToLog(FILE* fp, char* message);
FILE* openLogFile(char* filename);
int closeLogFile(FILE* fp);
int saveSynapseOutputFile(char* filename, void *obj, int duration, double dCpre, double dCpost, double dThetaD, double dThetaP, double dGammaD, double dGammaP, double dSigma, int iPreSpikeDelay, int iTau, int iTauC, double dRhoFixed, double poisson_param, long random_seed);


//int loadDataFile(char* filename, void *obj);
//int saveOutputFile(char* filename, void *obj);

#endif /*DATATOOLS_H_*/
