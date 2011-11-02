#include "GeneralIncludes.h"
#include "DataTools.h"
#include "LIFNeuron.h"
#include "Synapse.h"

#include "NumericalTools.h"
#include "SpikeTrains.h"
#include "ExternalCurrents.h"
#include "NetworkTopology.h"

signed int** network_memory_init(signed int **adj);
int assignNetworkTopology(signed int **adj);
int finalise(int status, LIFNeuron *lif, Synapse *syn, signed int **adj);

int main( int argc, char *argv[] )
{
    int i, j, t;
    char outfile[FILE_NAME_LENGTH];
    int syn_index;

    #ifndef _NO_DEBUG_
    printf("_NO_DEBUG_ not defined!!!\n");
    #else // _NO_DEBUG_
    printf("else option triggered, _NO_DEBUG_ must be defined\n");
    #endif // _NO_DEBUG_

    lif_time_of_last_save = -1;
    syn_time_of_last_save = -1;

    // Load LIFNeurons into memory
    LIFNeuron *lif = NULL;
    lif = neuron_checkpoint_init(argc, argv, lif);
    fflush(logfile);
    fflush(stdout);

    // Create adjacency matrix, default to -1
    signed int **adj = NULL;
    adj = network_memory_init(adj);

    // Assign network topology, no_synapses is defined by topology
    assignNetworkTopology(adj);

    // Debug network topology
    saveTopologyToFile(adj, outfilepattern);
    //printf("DEBUG: adjacency matrix\n");
    fflush(stdout);

    // Load Synapses into memory
    Synapse *syn = NULL;
    syn = synapse_checkpoint_init(argc, argv, syn);
    fflush(logfile);
    fflush(stdout);

    // Initialise checkpointing
    /* Checkpoint_init:
        1. Load parameters
        2. openLogFile() (now in load params)
        3. Reserve memory for array of synapses
        4. Load Reset values
        5. Reserve memory for members of Synapse(s)
    */

    // Load external driving currents
    loadExternalCurrent(lif);
    // Load pre- and post- synaptic spike times into arrays for each synapse
    loadInitialSpikeTimes(syn);

    // Setup raster output file
    createRasterOutputFileHeader(outfilepattern, siT, dCm, dRm, V_rest, V_reset, V_threshold, iRefracTime, initial_random_seed, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, current_transfer_delay, current_transfer_const, topology_id);

    // Main simulation loop
    fprintf(logfile, "Entering main simulation loop\n");
    printf("Entering main simulation loop\n");
    // Loop over discrete time steps up to simulation_duration
    for (t = siT; t < (simulation_duration-1); t++)
    {
        // Update each synapse
        for (i = 0; i < no_synapses; i++){
            printf("syn(%d) ", i);
            updateCalciumConcentration(&syn[i]);
            updateSynapticEfficacy(&syn[i]);
            printf("t: %d, c: %f, rho: %f\n", siT, syn[i].c[siT-syn_time_of_last_save], syn[i].rho[siT-syn_time_of_last_save]);
        }
        // Update each neuron
        for (i = 0; i < no_neurons; i++)
        {
            printf("lif(%d) ", i);
            updateNeuronMembraneVoltage(&lif[i]);
            printf("t: %d, V: %f, Spike: %d\n", siT, lif[i].V[siT-lif_time_of_last_save], lif[i].spikeT[siT-lif_time_of_last_save]);
        }
        // Transfer spikes from LIF.spikes to syn.pre/postT
        //TODO: share memory of lif.spikes and syn.pre/posT so this step isn't needed
        for (i = 0; i < no_neurons; i++){
            // Test each neuron to see if spike occurred
            if (lif[i].spikeT[siT + 1] == 1){
                // Across the post-synaptic neurons
                for (j = 0; j < no_neurons; j++){
                    // Transfer when spikes occur
                    if ( (syn_index = adj[i][j]) > -1 ){
                        //printf("DEBUG: spikes->preT transfer\n");
                        // Transfer spikes to corresponding Synapse preT array
                        syn[syn_index].preT[siT + 1] += 1;
                    }
                }
                // Across the pre-synaptic neurons
                for (j = 0; j < no_neurons; j++){
                    // Transfer spikes to synaptic arrays
                    if ( (syn_index = adj[j][i]) > -1){
                        //printf("DEBUG: spikes->postT transfer\n");
                        syn[syn_index].postT[siT + 1] += 1;
                    }
                }
            }
        }
        // Transfer currents
        if (siT > current_transfer_delay){
            for (i = 0; i < no_neurons; i++){
                // Test each neuron to see if spike occurred
                if (lif[i].spikeT[siT + 1 - current_transfer_delay] == 1){
                    // Output raster data
                    sprintf(outfile, outfilepattern, "raster", 0);
                    printf("writing (raster)...%s\n", outfile);
                    saveRasterOutput(outfile, i, (siT + 1 - current_transfer_delay));

                    //printf("DEBUG: Transfer spike as current! lif(%d)\n", i);
                    // Across the post-synaptic neurons
                    for (j = 0; j < no_neurons; j++){
                        // Transfer currents when spikes occur
                        if ( (syn_index = adj[i][j]) > -1 ){
                            // CONSIDER: is a linear current transfer fn reasonable?
                            lif[j].Isyn[siT + 1] += current_transfer_const * syn[syn_index].rho[siT + 1 - current_transfer_delay];
                        }
                    }
                }
            }
        }

        fflush(stdout);
        synapse_checkpoint_save(syn);
        neuron_checkpoint_save(lif);
        siT++;
    }
    printf("SIM OVER\n");
    synapse_checkpoint_save(syn);
    neuron_checkpoint_save(lif);
    for (i = 0; i < no_synapses; i++){
        printf("syn(%d) t: %d, c: %f, rho: %f\n", i, siT, syn[i].c[siT], syn[i].rho[siT]);
    }
    for (i = 0; i < no_neurons; i++)
    {
        printf("lif(%d) t: %d, V: %f, Spike: %d\n", i, siT, lif[i].V[siT], lif[i].spikeT[siT]);
    }
    fprintf(logfile, "Simulation complete\n");
    printf("Simulation complete\n");

    #ifndef _NO_DEBUG_
    // Debugging output after simulation has completed
    for (j = 0; j < (simulation_duration); j++){
        for (i = 0; i < no_synapses; i++){
            fprintf(logfile, "syn(%d).preT(%d): %u, postT(%d): %u, c: %f, rho: %f\n", i, j, syn[i].preT[j], j, syn[i].postT[j], syn[i].c[j], syn[i].rho[j]);
        }
    }
    for (j = 0; j < (simulation_duration); j++){
        for (i = 0; i < no_neurons; i++){
            fprintf(logfile, "lif(%d).Iext(%d): %f, Isyn(%d): %f, V: %f, Spike: %u\n", i, j, lif[i].Iext[j], j, lif[i].Isyn[j], lif[i].V[j], lif[i].spikeT[j]);
        }
    }
    fprintf(logfile, "siT: %d\n", siT);
    #endif // _NO_DEBUG_

    // Output to files loop
    if (!checkpointing)
    {
        for (i = 0; i < no_synapses; i++){
            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
            sprintf(outfile, outfilepattern, "syn", syn[i].ID);
            printf("writing...%s\n", outfile);
            saveSynapseOutputFile(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, initial_random_seed);
        }
        for (i = 0; i < no_neurons; i++){
            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
            sprintf(outfile, outfilepattern, "lif", lif[i].ID);
            printf("writing...%s\n", outfile);
            saveNeuronOutputFile(outfile, &lif[i], siT, dCm, dRm, V_rest, V_reset, V_threshold, iRefracTime, initial_random_seed);
        }
    }

    // Free memory and exit
    finalise(0, lif, syn, adj);
    return 0;
}


// Assign memory for adjacency matrix
signed int** network_memory_init(signed int **adj){
    int i, j;

    // Create adjacency matrix, default to 0
    adj = (signed int **) calloc( no_neurons, sizeof(signed int *) );

    if (adj == NULL){
        perror("Memory allocation failure (adj)\n");
        fprintf(logfile, "ERROR: Memory allocation failure (adj)\n");
    }
    else{
        for (i = 0; i < no_neurons; i++){
            adj[i] = calloc(no_neurons, sizeof(signed int));
            if (adj[i] == NULL){
                perror("Memory allocation failure (adj(row))\n");
                fprintf(logfile, "ERROR: Memory allocation failure (adj(row))\n");
            }
            else{
                for (j = 0; j < no_neurons; j++){
                    adj[i][j] = -1;
                }
            }
        }

        fprintf(logfile, "adj successfully assigned\n");
    }

    return adj;
}

// Create network topology
int assignNetworkTopology(signed int **adj){
    // Use topology functions defined in NetworkTopology.h
    no_synapses = topology_fn(adj, no_neurons, random_seed);

    //printf("DEBUG: no_synapses: %d\n", no_synapses);
    return 0;
}


// Free memory, close log file, ready to exit
int finalise(int status, LIFNeuron *lif, Synapse *syn, signed int **adj){
    int i;

    if (status == 0){
        synapse_finalise(0, syn);
        neuron_finalise(0, lif);
        for (i = 0; i < no_neurons; i++){
            free(adj[i]);
        }
        free(adj);
        closeLogFile(logfile);
        return 0;
    }
    else{
        fprintf(logfile, "An error occurred, exiting.");
        closeLogFile(logfile);
        return 1;
    }
}
