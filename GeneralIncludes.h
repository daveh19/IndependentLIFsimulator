#ifndef GENERAL_INCLUDES_H_
#define GENERAL_INCLUDES_H_

//#define DEBUG_MODE // not currently used
//#define DEBUG_MODE_MAIN // screen display of timestep based change of V and RHO
//#define DEBUG_MODE_NETWORK // print to file of connectivity statistics
//#define DEBUG_MODE_SPIKES // screen display of spike transfers
//#define DEBUG_MODE_SYNAPSE // screen display of synapse updates


#define FILE_NAME_LENGTH (50)
#define TEXT_BUFFER_LENGTH (100)
#define EPSILLON (0.0000001)

// Data reporters
#define RECORDER_NEURON_ID (201)
#define RECORDER_SYNAPSE_ID (10)

#define USE_GPU (1) /* 1=gpu, 0=cpu */
#define NETWORK_SEED (-14)
#define PARALLEL_SEED (0)
#define GAUSSIAN_SYNAPTIC_SEED (-12)
#define UNIFORM_SYNAPTIC_SEED (-10)

#define MAX_TIME_STEPS (300000) /*no of timesteps, each of size dt*/

// Network schema
#define NO_EXC (8000) /*(400)*/ /*(10000)*/
#define NO_INH (2000)
#define NO_LIFS (NO_EXC + NO_INH) 
#define CONNECTIVITY_PROBABILITY (0.05)

// Time step sizes and statistical bin widths
#define LIF_DT (0.0001) /* modify refrac time and calcium delay in tandem, also MAX_TIME_STEPS */
#define SYN_DT LIF_DT /*TODO: at a later stage I will have the synapse update more slowly than the lif*/
#define BIN_SIZE (0.1)

// Transfer voltages
#define J_EE (0.2) /*(0.1)*/
#define J_IE (0.2)
#define J_II (-0.8)
#define J_EI (-0.4)

#define RHO_FIXED (0.5) /*testing fixed transfer strength for Mean Field comparison*/

#define J_EXT (10.)
//TODO: separate drive to exc and inh pops

// LIF specific
#define LIF_V_INITIAL (-66.0)
#define LIF_V_REST (-70.0)
#define LIF_V_RESET (-68.0)
#define LIF_V_THRESHOLD (-54.0)
#define LIF_CM (0.001)
#define LIF_RM (20.0)
#define LIF_SIGMA (5)
#define LIF_REFRAC_TIME (200) /*timesteps*/

// Synapse model specific
#define SYN_RHO_INITIAL (1.0)
#define SYN_CA_INITIAL (0.0)
#define SYN_CALCIUM_DELAY (50) /*5*/ /*timesteps (needs to be modified when DT is modified above*/
#define SYN_GAMMA_P (725.085)
#define SYN_GAMMA_D (331.909)
#define SYN_THETA_P (1.3)
#define SYN_THETA_D (1.0)
#define SYN_SIGMA (3.35) /*(3.35)*/ /*3.35;*/ /*TODO: switch synapse noise back on*/
#define SYN_TAU (346.3615)
#define SYN_TAU_CA (0.0226936)
#define SYN_C_PRE (0.5617539)
#define SYN_C_POST (1.23964)



#include <stdio.h>
#include <errno.h> // Error numbers for file opening
//#include <limits.h> // System data limits
//#include <float.h> // Limits of system implementation of floating values
#include <stdlib.h> // malloc()
#include <fcntl.h> // manipulate file descriptors, opencl
//#include <strcmp.h>
#include <math.h> // supposedly for fmin()
#include <time.h> // to time main loop
#include <string.h> // strcpy() strcat()
//#include <ctype.h>
//#include <assert.h>
//#include <locale.h>
//#include <stddef.h>
#include <sys/types.h> // opencl
#include <sys/stat.h> // For mkdir()
//#include <sys/dir.h>

//from Mac OpenCL demo
#include <unistd.h>
#ifdef __APPLE__ //ifdef added from web suggestions
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#endif /*GENERAL_INCLUDES_H_*/
