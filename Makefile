#
# Simple makefile to compile the OpenCL cldemo.c
# on a Mac OSX with linking to the OpenCL library installed
#

CC = gcc
LD = gcc -lm

CFLAGS = -Wall -std=gnu99
CDEBUG =

LIBOCL = -L/Developer/SDKs/MacOSX10.6.sdk/System/Library/Frameworks/OpenCL.framework/Versions/A/Libraries
INCOCL = -I/System/Library/Frameworks/OpenCL.framework/Versions/A/Headers

INCCUDA = -I/usr/local/cuda/include/
LIBCUDA = -L/usr/local/cuda/ -L/usr/local/cuda/lib/ -l OpenCL

INCAMD = -I/opt/AMDAPP/include/
LIBAMD = -L/opt/AMDAPP/lib/x86_64/ -lOpenCL

SRCS = DataReporters.c HandleOpenCL.c NumericalTools.c main.c
OBJS = *.o
 
#EXE = kernel/network_sim
EXE = build/Debug/net_sim_2

midway:
				$(LD) $(CFLAGS) -lOpenCL $(SRCS) -o $(EXE)
                
lynch:
				$(LD) $(CFLAGS) $(INCCUDA) $(LIBCUDA) $(SRCS) -o $(EXE)
                
bsd:
				$(LD) $(CFLAGS) $(INCAMD) $(LIBAMD) $(SRCS) -o $(EXE)

all: $(EXE)

$(OBJS): $(SRCS)
				$(CC) $(CFLAGS) $(INCOCL) -I/usr/include -c $(SRCS)

$(EXE): $(OBJS)
				$(LD) -L/usr/local/lib $(OBJS) $(LIBOCL) -o $(EXE) -framework OpenCL

clean:
				rm -f $(OBJS) *~
				#clear
