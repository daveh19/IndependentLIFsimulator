#include "GeneralIncludes.h"
#include "Synapse.h"
#include "DataTools.h"


#include "NumericalTools.h"
#include "SpikeTrains.h"

int main( int argc, char *argv[] ){
    int i, j, t;
    char outfile[FILE_NAME_LENGTH];

    // Initialise checkpointing
/* Checkpoint_init:
    1. Load parameters
    2. openLogFile() (now in load params)
    3. Reserve memory for array of synapses
    4. Load Reset values
    5. Reserve memory for members of Synapse(s)
*/
    Synapse *syn;
    syn = checkpoint_init(argc, argv, syn);
    fflush(logfile);

//    for (k = 0; i < no_synapses; i++){
//        fprintf(logfile, "DEBUG:: main: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
//    }
//    printf("DEBUG:: MEM TEST\n");
//    printf("syn.ID: %d\n", syn[0].ID);
    fflush(stdout);



    // Load pre- and post- synaptic spike times into arrays for each synapse
    loadInitialSpikeTimes(syn);

    // Main simulation loop
    fprintf(logfile, "Entering main simulation loop\n");
    printf("Entering main simulation loop\n");
    // Loop over discrete time steps up to simulation_duration
    for (t = siT; t < (simulation_duration-1); t++){
        //checkpoint_save(syn);
        // Update each synapse
        for (i = 0; i < no_synapses; i++){
            printf("syn(%d) ", i);
            updateCalciumConcentration(&syn[i]);
            updateSynapticEfficacy(&syn[i]);
            printf("t: %d, c: %f, rho: %f\n", siT, syn[i].c[siT-time_of_last_save], syn[i].rho[siT]);
        }
        checkpoint_save(syn);
        siT++;
    }
    printf("DEBUG:: SIM OVER\n");
    checkpoint_save(syn);
    for (i = 0; i < no_synapses; i++){
        printf("syn(%d) t: %d, c: %f, rho: %f\n", i, siT, syn[i].c[siT], syn[i].rho[siT]);
    }
    fprintf(logfile, "Simulation complete\n");
    printf("Simulation complete\n");

    // Debugging output after simulation has completed
    for (j = 0; j < (simulation_duration); j++){
        for (i = 0; i < no_synapses; i++){
            fprintf(logfile, "syn(%d).preT(%d): %u, postT(%d): %u, c: %f, rho: %f\n", i, j, syn[i].preT[j], j, syn[i].postT[j], syn[i].c[j], syn[i].rho[j]);
        }
    }
    fprintf(logfile, "siT: %d\n", siT);

    // Output to files loop
    if (!checkpointing){
        for (i = 0; i < no_synapses; i++){
            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
            sprintf(outfile, outfilepattern, syn[i].ID);
            printf("writing...%s\n", outfile);
            saveSynapseOutputFile(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, poisson_param, initial_random_seed);
        }
    }

    // Free memory and exit
    return finalise(0, syn);
}


// Calculate synaptic efficacy for next time step
void updateSynapticEfficacy(Synapse *syn){
    double rho, drho, minTheta, rand_no, noise;
    rho = (*syn).rho[siT];
    drho = (-rho * (1.0 - rho) * (dRhoFixed - rho)) + (dGammaP * (1 - rho) * h(syn, dThetaP)) - (dGammaD * rho * h(syn, dThetaD));

    // Add noise
    minTheta = fmin(dThetaP, dThetaD);
    if (h(syn, minTheta) > 0){ // Noise is on
        rand_no = (double) gasdev(&random_seed);
        noise = dSigma * sqrt(iTau) * rand_no;
        printf("\nNoise is active, rand_no: %f, noise: %f\n", rand_no, noise);
        drho += noise;
    }

    drho /= (double)iTau;
    if ( (rho + drho) > 0 ){
        (*syn).rho[siT + 1] = rho + drho; // Euler forward method
    }
    else{
        (*syn).rho[siT + 1] = 0.0;
    }
}


// Simple Heaviside implemenation for comparing calcium
// concentration with a threshold value
BOOL h(Synapse *syn, double theta){
    if ( (*syn).c[siT] < theta)
        return 0;
    else
        return 1;
}

// Calculate synaptic calcium concentration for next time step
void updateCalciumConcentration(Synapse *syn){
    double c, dc;
    c = (*syn).c[siT];
    dc = (-c / (double)iTauC) + calciumFromPreSynapticSpikes(syn) + calciumFromPostSynapticSpikes(syn);
    (*syn).c[siT + 1] = c + dc; // Euler forward method
}


// Calculate contribution to next synaptic calcium concentration
// from pre-synaptic spikes
// Note: there is a delay iPreSpikeDelay before calcium from a
// pre-synaptic spike enters the synaptic cleft
double calciumFromPreSynapticSpikes(Synapse *syn){
    double d;

    printf("preT: %u ", (*syn).preT[siT]);

    if (siT < iPreSpikeDelay){
        d = 0.0;
    }
    else if( (siT >= iPreSpikeDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iPreSpikeDelay]) * dCpre;
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in calciumFromPreSynapticSpikes()");
    }

    return d;
}


// Calculate contribution to next synaptic calcium concentration
// from post-synaptic spikes
double calciumFromPostSynapticSpikes(Synapse *syn){
    double d;
    printf("postT: %u ", (*syn).postT[siT]);
    d = ((double) (*syn).postT[siT]) * dCpost;
    return d;
}


// Setup spike times (hard-coded version)
// preT[i] = 1 means a spike occurs at time i
// preT[i] = 0 implies no spike at time i
void loadInitialSpikeTimes(Synapse *syn){
    int i;
    fprintf(logfile, "Initialising spike times\n");
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
    for (i = 0; i < no_synapses; i++){
        (*train_fn)(syn[i].preT, syn[i].postT, simulation_duration);
    }
    fprintf(logfile, "Spike times initialised\n");
    //fflush(logfile);
}


void synapse_memory_init(Synapse *syn){
    int i;
    double * local_c;
    double * local_rho;
    unsigned int * local_preT;
    unsigned int * local_postT;
    //Synapse * local_synapse;
    fprintf(logfile, "Synapse simulator initialising.\n");

    for (i = 0; i < no_synapses; i++){
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
        (syn[i]).ID = siID;
        siID++;
        fprintf(logfile, "Set synaptic id to: %d\n", (syn[i]).ID);

        // Memory allocation for c(t) array
        local_c = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_c == NULL){
            perror("Memory allocation failure (c)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (c)\n");
        }
        else{//removed (*syn) to allow for array based syn[0]
            (syn[i]).c = local_c;
            syn[i].c[0] = initial_c; // TODO: check if this hardcoded 0 is ok
            fprintf(logfile, "syn(%d).c successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
        }
        // Memory allocation for rho(t) array
        local_rho = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_rho == NULL){
            perror("Memory allocation failure (rho)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (rho)\n");
        }
        else{
            (syn[i]).rho = local_rho;
            //syn[i].rho[0] = initial_rho;
            fprintf(logfile, "syn(%d).rho successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).rho(0): %lf\n", i, syn[i].rho[0]);
        }

        // Memory allocation for preT(t) array
        // CONSIDER: using calloc instead of malloc for spike time arrays (defaults to 0)
        local_preT = (unsigned int *) malloc( (simulation_duration) * sizeof(unsigned int) );
        if (local_preT == NULL){
            perror("Memory allocation failure (preT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (preT)\n");
        }
        else{
            (syn[i]).preT = local_preT;
            fprintf(logfile, "syn(%d).preT successfully assigned\n", i);
            //(syn[i]).preT[0] = 99;  //
            //fprintf(logfile, "DEBUG:: syn(%d).preT[0] is %d\n", i, syn[i].preT[0]);
            //fflush(logfile);
        }
        // Memory allocation for postT(t) array
        local_postT = (unsigned int *) malloc( (simulation_duration) * sizeof(unsigned int) );
        if (local_postT == NULL){
            perror("Memory allocation failure (postT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (postT)\n");
        }
        else{
            (syn[i]).postT = local_postT;
            fprintf(logfile, "syn(%d).postT successfully assigned\n", i);
        }
    }
    fprintf(logfile, "Initialisation of simulator complete\n");
}


int finalise(int status, Synapse *syn){
    int i;
    if (status == 0){
        fprintf(logfile, "Synapse simulator exiting successfully\n");
        for (i = 0; i < no_synapses; i++){
            free((syn[i]).c);
            free((syn[i]).rho);
            free((syn[i]).preT);
            free((syn[i]).postT);
        }
        free(syn);
        fprintf(logfile, "Memory freed\n");
        fprintf(logfile, "Exiting\n");
        closeLogFile(logfile);
        return 0;
    }
    else{
        fprintf(logfile, "An error occurred: exiting\n");
        closeLogFile(logfile);
        return 1;
    }
}
