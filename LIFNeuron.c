#include "GeneralIncludes.h"
#include "LIFNeuron.h"
#include "DataTools.h"


#include "NumericalTools.h"
#include "SpikeTrains.h"

//int main( int argc, char *argv[] ){
//    int i, j, t;
//    char outfile[FILE_NAME_LENGTH];
//
//    // Initialise checkpointing
///* Checkpoint_init:
//    1. Load parameters
//    2. openLogFile() (now in load params)
//    3. Reserve memory for array of synapses
//    4. Load Reset values
//    5. Reserve memory for members of Synapse(s)
//*/
//    LIFNeuron *lif;
//    lif = neuron_checkpoint_init(argc, argv, lif);
//    fflush(logfile);
//
////    for (k = 0; i < no_synapses; i++){
////        fprintf(logfile, "DEBUG:: main: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
////    }
////    printf("DEBUG:: MEM TEST\n");
////    printf("syn.ID: %d\n", syn[0].ID);
//    fflush(stdout);
//
//
//
//    // Load pre- and post- synaptic spike times into arrays for each synapse
//    loadExternalCurrent(lif);
//
//    // Main simulation loop
//    fprintf(logfile, "Entering main simulation loop\n");
//    printf("Entering main simulation loop\n");
//    // Loop over discrete time steps up to simulation_duration
//    for (t = siT; t < (simulation_duration-1); t++){
//        //checkpoint_save(syn);
//        // Update each neuron
//        for (i = 0; i < no_neurons; i++){
//            printf("lif(%d) ", i);
//            //updateCalciumConcentration(&syn[i]);
//            //updateSynapticEfficacy(&syn[i]);
//            updateNeuronMembraneVoltage(&lif[i]);
//            printf("t: %d, V: %f, Spike: %d\n", siT, lif[i].V[siT-lif_time_of_last_save], lif[i].spikeT[siT]);
//        }
//        neuron_checkpoint_save(lif);
//        siT++;
//    }
//    printf("DEBUG:: SIM OVER\n");
//    neuron_checkpoint_save(lif);
//    for (i = 0; i < no_neurons; i++){
//        printf("lif(%d) t: %d, V: %f, Spike: %d\n", i, siT, lif[i].V[siT], lif[i].spikeT[siT]);
//    }
//    fprintf(logfile, "Simulation complete\n");
//    printf("Simulation complete\n");
//
//    // Debugging output after simulation has completed
//    for (j = 0; j < (simulation_duration); j++){
//        for (i = 0; i < no_neurons; i++){
//            fprintf(logfile, "lif(%d).Iext(%d): %f, Isyn(%d): %f, V: %f, Spike: %u\n", i, j, lif[i].Iext[j], j, lif[i].Isyn[j], lif[i].V[j], lif[i].spikeT[j]);
//        }
//    }
//    fprintf(logfile, "siT: %d\n", siT);
//
//    // Output to files loop
//    if (!checkpointing){
//        for (i = 0; i < no_neurons; i++){
//            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
//            sprintf(outfile, neuron_outfilepattern, lif[i].ID);
//            printf("writing...%s\n", outfile);
//            saveNeuronOutputFile(outfile, &lif[i], siT, dCm, dRm, V_rest, V_reset, V_threshold, iRefracTime, initial_random_seed);
//        }
//    }
//
//    // Free memory and exit
//    return neuron_finalise(0, lif);
//}


// Calculate neuronal membrane dynamics for next time step
void updateNeuronMembraneVoltage(LIFNeuron *lif){
    double v, dv, newV, noise;
    int s;
    int i;

    double tau_m = dCm * dRm;
    float gauss;

    v = (*lif).V[siT];
    dv = 0;
    s = 0;
    gauss = gasdev(&random_seed);

    for(i = fmax((siT-iRefracTime+1),0); i < siT; i++){
        // Count number of spikes during refractory period, but preceeding current time point
        s += (*lif).spikeT[i];
    }
    if(s > 0){
        // During refractory period, but a spike has not 'just' occurred
        // V should have already been reset on previous time step
        v = (*lif).V[siT];
        // Apply leak current
        dv = (-(v - V_rest)/(dCm * dRm));
        // Apply noise
        noise = sqrt(dDt/tau_m) * lifSigma * gauss;
        newV = v + (dv * dDt) + noise; //CONSIDER: should I be applying noise here?
        // Apply lower threshold on Vm
        if (newV < V_rest){
            newV = V_rest;
        }
        // No spike occurred, append that to the list
        (*lif).spikeT[siT + 1] = 0;
        #ifndef _NO_DEBUG_
        fprintf(logfile, "ref--t: %d dv: %f newV: %f\n", siT, dv, newV);
        #endif // _NO_DEBUG_
    }
    else if ((*lif).spikeT[siT] > 0){
        // A spike has just occurred
        // During refractory period, so reset Vm to Vreset then apply leak
        if (iRefracTime > 0){
            v = V_reset;
            // Apply leak current
            dv = (-(v - V_rest)/(dCm * dRm));
            // Apply noise
            noise = sqrt(dDt/tau_m) * lifSigma * gauss;
            newV = v + (dv * dDt) + noise; //CONSIDER: should I be applying noise here?
            // Apply lower threshold on Vm
            if (newV < V_rest){
                newV = V_rest;
            }
            // No spike occrred, append that to the list
            (*lif).spikeT[siT + 1] = 0;
            #ifndef _NO_DEBUG_
            fprintf(logfile, "refrac time, dv: %f newV: %f\n", dv, newV);
            #endif // _NO_DEBUG_
        }
        // A spike has just occurred, but refractory period is 0 so reset Vm to Vreset then immediately calculate and add dv
        else{
            v = V_reset;
            // Apply leak current
            dv = (-(v - V_rest)/(dCm * dRm));
            // If external input current exists then apply it to Vm
            dv += ((*lif).Iext[siT] / dCm);
            // If synaptic input current exists apply it to Vm
            dv += ((*lif).Isyn[siT] / dCm);
            // Apply noise
            noise = sqrt(dDt/tau_m) * lifSigma * gauss;
            newV = v + (dv * dDt) + noise;
            // Apply lower threshold on Vm
            if (newV < V_rest){
                newV = V_rest;
            }
            if (newV > V_threshold){
                #ifndef _NO_DEBUG_
                fprintf(logfile, "(b) Spike occurred, t: %u newV: %f", siT, newV);
                #endif // _NO_DEBUG_
                (*lif).spikeT[siT + 1] = 1;
                newV = -20; // Make the spike nice and obvious on plots #CONSIDER: find a Vm value for spikes
            }
            else{
                (*lif).spikeT[siT + 1] = 0;
            }
            #ifndef _NO_DEBUG_
            fprintf(logfile, "(b) Neuron(%d) dv: %f newV: %f\n", (*lif).ID, dv, newV);
            #endif // _NO_DEBUG_
        }
    }
    else{
        // Not in refractory period, no spike has just occurred, so update V according to dV
        v = (*lif).V[siT];
        // Apply leak current
        dv = (-(v - V_rest) / (dCm * dRm));
        // Apply external input current
        dv += ((*lif).Iext[siT] / dCm);
        // Apply synaptic input current
        dv += ((*lif).Isyn[siT] / dCm);
        // Apply noise
        noise = sqrt(dDt/tau_m) * lifSigma * gauss;
        newV = v + (dv * dDt) + noise;
        // Apply lower threshold to Vm
        if (newV < V_rest){
            newV = V_rest;
        }
        if (newV > V_threshold){
            #ifndef _NO_DEBUG_
            fprintf(logfile, "(a) Spike occurred, t: %u newV: %f", siT, newV);
            #endif // _NO_DEBUG_
            (*lif).spikeT[siT + 1] = 1;
            newV = -20; // Make the spike nice and obvious on plots #CONSIDER: find a Vm value for spikes
        }
        else{
            (*lif).spikeT[siT + 1] = 0;
        }
        #ifndef _NO_DEBUG_
        fprintf(logfile, "(a) Neuron(%d) dv: %f newV: %f\n", (*lif).ID, dv, newV);
        #endif // _NO_DEBUG_
    }
    (*lif).V[siT + 1] = newV;
}


//// Simple Heaviside implemenation for comparing calcium
//// concentration with a threshold value
//BOOL h(LIFNeuron *lif, double theta){
//    if ( (*lif).V[siT] < theta)
//        return 0;
//    else
//        return 1;
//}


void loadExternalCurrent(LIFNeuron *lif){
    int i;
    fprintf(logfile, "Initialising external currents\n");
//    fflush(logfile);
//    fprintf(logfile, "DEBUG:: syn(%d).preT[0] is %d\n", 0, (*syn).preT[0]);
//    fflush(logfile);
//    syn[0].preT[0] = 1;
//    syn[0].postT[0] = 0;
//    printf("DEBUG:: first spikes\n");
//    for (i = 1; i < simulation_duration; i++){
//        syn[0].preT[i] = 0;
//        syn[0].postT[i] = 0;
//    }
    for (i = 0; i < no_neurons; i++){
        (*current_fn)(lif[i].Iext, simulation_duration);
    }
    fprintf(logfile, "External currents initialised\n");
    //fflush(logfile);
}


void lif_memory_init(LIFNeuron *lif){
    int i;
    double * local_V;
    double * local_Iext;
    double * local_Isyn;
    unsigned int * local_spikeT;
    //Synapse * local_synapse;
    fprintf(logfile, "Neuron simulator initialising.\n");

    for (i = 0; i < no_neurons; i++){
//        // Memory allocation for each synapse
//        local_synapse = (Synapse *) malloc( sizeof(Synapse) );
//        if (local_synapse == NULL){
//            perror("Memory allocation error (Synapse)\n");
//            fprintf(logfile, "ERROR: Memory allocation failure (Synapse)\n");
//        }
//        else{
//            (syn[i]) = local_synapse;
//            fprintf(logfile, "syn(%d) successfully assigned\n", i);
//        }
        // Set synapse ID
        (lif[i]).ID = siID;
        siID++;
        fprintf(logfile, "Set neuron id to: %d\n", (lif[i]).ID);

        // Memory allocation for V(t) array
        local_V = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_V == NULL){
            perror("Memory allocation failure (V)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (V)\n");
        }
        else{//removed (*syn) to allow for array based syn[0]
            (lif[i]).V = local_V;
            lif[i].V[0] = initial_v; // TODO: check if this hardcoded 0 is ok
            fprintf(logfile, "lif(%d).V successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
        }
        // Memory allocation for Iext(t) array
        local_Iext = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_Iext == NULL){
            perror("Memory allocation failure (Iext)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (Iext)\n");
        }
        else{
            (lif[i]).Iext = local_Iext;
            //syn[i].rho[0] = initial_rho;
            fprintf(logfile, "lif(%d).Iext successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).rho(0): %lf\n", i, syn[i].rho[0]);
        }
        // Memory allocation for Isyn(t) array
        local_Isyn = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_Isyn == NULL){
            perror("Memory allocation failure (Isyn)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (Isyn)\n");
        }
        else{
            (lif[i]).Isyn = local_Isyn;
            //syn[i].rho[0] = initial_rho;
            fprintf(logfile, "lif(%d).Isyn successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).rho(0): %lf\n", i, syn[i].rho[0]);
        }

        // Memory allocation for spikeT(t) array
        // CONSIDER: using calloc instead of malloc for spike time arrays (defaults to 0)
        local_spikeT = (unsigned int *) malloc( (simulation_duration) * sizeof(unsigned int) );
        if (local_spikeT == NULL){
            perror("Memory allocation failure (spikeT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (spikeT)\n");
        }
        else{
            (lif[i]).spikeT = local_spikeT;
            fprintf(logfile, "lif(%d).spikeT successfully assigned\n", i);
            //(syn[i]).preT[0] = 99;  //
            //fprintf(logfile, "DEBUG:: syn(%d).preT[0] is %d\n", i, syn[i].preT[0]);
            //fflush(logfile);
        }
    }
    fprintf(logfile, "Initialisation of simulator complete\n");
}


int neuron_finalise(int status, LIFNeuron *lif){
    int i;
    if (status == 0){
        fprintf(logfile, "Neuron simulator exiting successfully\n");
        for (i = 0; i < no_neurons; i++){
            free((lif[i]).V);
            free((lif[i]).Iext);
            free((lif[i]).Isyn);
            free((lif[i]).spikeT);
        }
        free(lif);
        fprintf(logfile, "LIFNeuron memory freed\n");
        fprintf(logfile, "Exiting\n");
        //closeLogFile(logfile);
        return 0;
    }
    else{
        fprintf(logfile, "An error occurred: exiting\n");
        //closeLogFile(logfile);
        return 1;
    }
}
