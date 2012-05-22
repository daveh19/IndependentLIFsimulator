#ifndef GENERAL_INCLUDES_H_
#define GENERAL_INCLUDES_H_

#define FILE_NAME_LENGTH (50)
#define TEXT_BUFFER_LENGTH (100)
#define EPSILLON (0.0000001)

#define RECORDER_NEURON_ID (100)
#define RECORDER_SYNAPSE_ID (0)

#define USE_GPU (1) /* 1=gpu, 0=cpu */
#define NETWORK_SEED (-14)

#define MAX_TIME_STEPS (3000)

#define NO_EXC (8000) /*(400)*/ /*(10000)*/ /*(400)*/
#define NO_INH (2000)
#define NO_LIFS (NO_EXC + NO_INH) 
#define CONNECTIVITY_PROBABILITY (0.05)

#define LIF_DT (0.001)
#define BIN_SIZE (0.1)

#define CALCIUM_DELAY (5)

#define J_EE (0.1)
#define J_IE (0.2)
#define J_II (-0.8)
#define J_EI (-0.4)

#define J_EXT (10.)
//TODO: separate drive to exc and inh pops


#include <stdio.h>
#include <errno.h> // Error numbers for file opening
//#include <limits.h> // System data limits
//#include <float.h> // Limits of system implementation of floating values
#include <stdlib.h> // malloc()
#include <fcntl.h> // manipulate file descriptors, opencl
//#include <strcmp.h>
#include <math.h> // supposedly for fmin()
//#include <time.h>
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
