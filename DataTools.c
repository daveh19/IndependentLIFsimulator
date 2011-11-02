#include "GeneralIncludes.h"
#include "DataTools.h"
//#include "Synapse.h"
#include "SpikeTrains.h"


/* Checkpoint_init:
    1. Load parameters
    2. openLogFile() (now in load params())
    3. Reserve memory for array of synapses
    4. Load Reset values
    5. Reserve memory for members of Synapse(s)
*/
Synapse* checkpoint_init(int argc, char *argv[], Synapse *syn){
    char *checkpoint_filename, *restart;//, old_checkpoint_filename[100];
    FILE *checkpoint_fp;
    char outfile[FILE_NAME_LENGTH];
    int i;

    // Setup global variable of whether checkpointing is enabled
    if (! (checkpoint_filename = getenv("SYNAPSE_CHECKFILE")) ){
        checkpointing = 0;
    }
    else{
        checkpointing = 1;
    }

    if (! (restart = getenv("SYNAPSE_RESTART")) ){
        printf("No previous checkpoint set, starting afresh.\n");
        loadSimulationParameters(argc, argv);
        // Memory allocation for array of synapses
        syn = (Synapse *) malloc( no_synapses * sizeof(Synapse));
        if (syn == NULL){
            perror("Memory allocation failure (syn array)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (syn array)\n");
        }
        else{
            fprintf(logfile, "syn array successfully assigned %d synapses\n", no_synapses);
        }
//        fflush(logfile);
        synapse_memory_init(syn);
        for (i = 0; i < no_synapses; i++){
            syn[i].rho[siT] = initial_rho;
            syn[i].c[siT] = initial_c;
        }
//        for (i = 0; i < no_synapses; i++){
//            fprintf(logfile, "DEBUG:: checkpoint_init(1): syn(%d).c(0): %lf\n", i, syn[i].c[0]);
//        }

        // For checkpointing, must add headers separately to output files
        if ( checkpointing ){
            for (i = 0; i < no_synapses; i++){
                //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
                sprintf(outfile, outfilepattern, syn[i].ID);
                printf("writing...%s\n", outfile);
                createOutputFileHeader(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, poisson_param, initial_random_seed);
            }
        }
    }
    else{
        if (! (checkpoint_filename = getenv("SYNAPSE_CHECKFILE")) ){
            fprintf(stderr, "Error: no checkpoint file for the restart job.\n");
            exit(1);
        }
        else{
            printf("Restarting the job from %s.\n", checkpoint_filename);
            if(! (checkpoint_fp = fopen(checkpoint_filename, "r")) ){
                perror(checkpoint_filename);
                exit(2);
            }
            else{
                loadSimulationParameters(argc, argv);
                // Memory allocation for array of synapses
                syn = (Synapse *) malloc( no_synapses * sizeof(Synapse));
                if (syn == NULL){
                    perror("Memory allocation failure (syn array)\n");
                    fprintf(logfile, "ERROR: Memory allocation failure (syn array)\n");
                }
                else{
                    fprintf(logfile, "syn array successfully assigned\n");
                }
                synapse_memory_init(syn);

                // Attempt to load from checkpoint file
                if ( (checkpoint_load(checkpoint_fp, syn) > 0) ){
                    fprintf(stderr, "Error reading from restart file %s.\n", checkpoint_filename);
                    exit(3);
                }
                else{
                    printf("Reading from restart file succeeded.\n");
                }
                fclose(checkpoint_fp);
            }
        }
    }

    //time_of_last_save = siT;

//    for (i = 0; i < no_synapses; i++){
//        for (j = 0; j < simulation_duration; j++)
//            fprintf(logfile, "DEBUG:: checkpoint_init(2): syn(%d).c(%d): %lf\n", i, j, syn[i].c[j]);
//    }
//    fflush(logfile);

    return syn;
}


// On failure return 0, otherwise return (meaningful?) value > 0
// call loadSimulationParameters() and synapse_memory_init() before checkpoint_load()
int checkpoint_load(FILE *checkpoint_fp, Synapse *syn){
    char paramName[FILE_NAME_LENGTH];
    double paramValue;
    int syn_index;
    int i;

    printf("Loading from checkpoint file.\n");
    while (!feof(checkpoint_fp)) {
        if (fscanf(checkpoint_fp, "%s %lf\n", paramName, &paramValue) != 2){
            fprintf(logfile, "Error reading previous spike times (a)\n");
            break;
        }
        printf("Read(a): %s %f\n", paramName, paramValue);
        if (!strcmp(paramName, "SIT")){
            printf("Previous siT: %d\n", (int) paramValue);
            siT = (int) paramValue;
            time_of_last_save = siT;
        }
        else if (!strcmp(paramName, "RANDOM_SEED")){
            printf("Previous random_seed: %ld\n", (long) paramValue);
            random_seed = (int) paramValue;
        }
        else if (!strcmp(paramName, "SYNAPSE")){
            syn_index = (int) paramValue;
            printf("Loading Synapse(%d)\n", syn_index);
            for (i = 0; i < 2; i++){
                if (fscanf(checkpoint_fp, "%s %lf\n", paramName, &paramValue) != 2){
                    fprintf(logfile, "Error reading previous spike times (b)\n");
                    break;
                }
                printf("Read(b): %s %f\n", paramName, paramValue);
                if (!strcmp(paramName, "LAST_C")){
                    printf("Previous c(t): %f\n", paramValue);
                    initial_c = paramValue;
                    syn[syn_index].c[siT] = initial_c;
    //                printf("DEBUG:: syn[syn_index].c[siT]: %lf\n", syn[syn_index].c[siT]);
                }
                else if (!strcmp(paramName, "LAST_RHO")){
                    printf("Previous rho(t): %f\n", paramValue);
                    initial_rho = paramValue;
                    syn[syn_index].rho[siT] = initial_rho;
                }
            }
            // Attempt to read spike times
            if (fscanf(checkpoint_fp, "%s ", paramName) != 1){
                fprintf(logfile, "Error reading previous spike times (c)\n");
                break;
            }
            printf("Read(d): %s\n", paramName);
            if (!strcmp(paramName, "SPIKES")){
                for (i = 0; i < (iPreSpikeDelay+1); i++){
                    if (fscanf(checkpoint_fp, "%lf ", &paramValue) != 1){
                        fprintf(logfile, "Error reading previous spike times (d)\n");
                        break;
                    }
                    printf("Read(e): %d\n", (int) paramValue);
                    syn[syn_index].preT[(siT-iPreSpikeDelay+i)] = (int) paramValue;
                    //printf("DEBUG:: syn[syn_index].preT[(siT-iPreSpikeDelay+i)]: %d\n", syn[syn_index].preT[(siT-iPreSpikeDelay+i)]);
                }
                fscanf(checkpoint_fp, "\n");
            }
//            if (fscanf(checkpoint_fp, "%s %lf\n", paramName, &paramValue) != 2)
//                break;
//            printf("Read(c): %s %f\n", paramName, paramValue);
//            if (!strcmp(paramName, "LAST_C")){
//                printf("Previous c(t): %f\n", paramValue);
//                initial_c = paramValue;
//                syn[syn_index].c[siT] = initial_c;
////                printf("DEBUG:: syn[syn_index].c[siT]: %lf\n", syn[syn_index].c[siT]);
//            }
//            else if (!strcmp(paramName, "LAST_RHO")){
//                printf("Previous rho(t): %f\n", paramValue);
//                initial_rho = paramValue;
//                syn[syn_index].rho[siT] = initial_rho;
////                printf("DEBUG:: syn[syn_index].rho[siT]: %lf\n", syn[syn_index].rho[siT]);
//            }
        }
    }
    //fclose(checkpoint_fp);
    printf("DEBUG:: exiting checkpoint_load()\n");
//    for(i = 0; i < no_synapses; i++){
//        for (j = 0; j < simulation_duration; j++){
//            printf("DEBUG:: syn(%d).preT(%d) %d\n", i, j, syn[i].preT[j] );
//        }
//    }
    return 0;
}


int checkpoint_save(Synapse *syn){
    char *checkpoint_name, *restart, old_checkpoint_name[FILE_NAME_LENGTH];
    FILE *checkpoint_fp;
    char outfile[FILE_NAME_LENGTH];
    int i, j;

    if (! (checkpoint_name = getenv ("SYNAPSE_CHECKFILE"))) {
        fprintf (logfile, "Checkpointing not enabled...skipping\n");
    }
    else {
        // Rename old checkpoint file
        if ( (restart = getenv("SYNAPSE_RESTART")) ) {
            strcpy (old_checkpoint_name, checkpoint_name);
            strcat (old_checkpoint_name, ".old");
            fprintf (logfile, "Renaming old checkpoint restart file to %s\n", old_checkpoint_name);
            if (0 > rename (checkpoint_name, old_checkpoint_name)) {
                perror (old_checkpoint_name);
                exit (4);
            }
        }
        // Save progress so far to output file
        fprintf(logfile, "Saving progress so far to %s\n", outfilearray);
        // Output to files loop
        for (i = 0; i < no_synapses; i++){
            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
            sprintf(outfile, outfilepattern, syn[i].ID);
            printf("writing...%s\n", outfile);
            saveSynapseProgressToFile(outfile, &syn[i], siT);
        }
        // Save restart data to checkpoint file
        fprintf (logfile, "Saving checkpoint data on %s\n", checkpoint_name);
        if (! (checkpoint_fp = fopen (checkpoint_name, "w"))) {
            perror (checkpoint_name);
            exit (5);
        }
        else {
            fprintf(checkpoint_fp, "SIT %d\n", siT);
            fprintf(checkpoint_fp, "RANDOM_SEED %ld\n", random_seed);
            for (i = 0; i < no_synapses; i++){
                fprintf(checkpoint_fp, "SYNAPSE %d\n", i);
                fprintf(checkpoint_fp, "LAST_C %lf\n", syn[i].c[siT]);
                fprintf(checkpoint_fp, "LAST_RHO %lf\n", syn[i].rho[siT]);
                fprintf(checkpoint_fp, "SPIKES ");
                for (j = (siT-iPreSpikeDelay); j <= siT; j++){
                    if (j >= 0){
                        fprintf(checkpoint_fp, " %d", syn[i].preT[j]);
                    }
                }
                fprintf(checkpoint_fp, "\n");
            }
            fclose (checkpoint_fp);
        }
    }

    return 0;
}


int saveSynapseProgressToFile(char* filename, void *obj, int end_time ){
    FILE *fp;
    int i;

    Synapse * syn = (Synapse *)obj;

    fprintf(logfile, "Opening file for saving: %s\n", filename);
    fp = fopen(filename, "a");
    if (fp == NULL){
        perror("Error in saveSynapseProgressToFile()");
    }
    else{
        fprintf(logfile, "Saving to file: %s\n", filename);
        printf("DEBUG:: before\n");
        for (i = (time_of_last_save+1); i <= end_time; i++){
            fprintf(fp, "%d %f %f %u %u\n", i, (*syn).rho[i], (*syn).c[i], (*syn).preT[i], (*syn).postT[i]);
            printf("DEBUG: saving progress %d %f %f %u %u\n", i, (*syn).rho[i], (*syn).c[i], (*syn).preT[i], (*syn).postT[i]);
        }
        //fprintf(fp,"\n\n");
        printf("DEBUG:: after\n");
        fclose(fp);
        fprintf(logfile, "Completed saving\n");
        time_of_last_save = end_time;
    }

    return 0;
}


int createOutputFileHeader(char* filename, void *obj, int duration, double dCpre, double dCpost, double dThetaD, double dThetaP, double dGammaD, double dGammaP, double dSigma, int iPreSpikeDelay, int iTau, int iTauC, double dRhoFixed, double poisson_param, long initial_random_seed){
    FILE *fp;

    Synapse * syn = (Synapse *)obj;

    fprintf(logfile, "Opening file for saving: %s\n", filename);
    fp = fopen(filename, "a");
    if (fp == NULL){
        perror("Error in saveSynapseOutputFile()");
    }
    else{
        fprintf(logfile, "Saving to file: %s\n", filename);

        fprintf(fp, "\n\n# Params:\n#    Cpre: %f, Cpost: %f, thetaD: %f, thetaP: %f, gammaD: %f, gammaP: %f, sigma: %f, Delay: %d, tau: %d, tauC: %d, rhoF: %f, poisson param: %f, seed: %ld\n", dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, poisson_param, initial_random_seed);
        fprintf(fp,"#\n# Synaptic output, Synapse(%d):\n# t rho c preT postT\n", (*syn).ID);

        fclose(fp);
        fprintf(logfile, "Completed saving\n");
    }

    return 0;
}


int saveSynapseOutputFile(char* filename, void *obj, int duration, double dCpre, double dCpost, double dThetaD, double dThetaP, double dGammaD, double dGammaP, double dSigma, int iPreSpikeDelay, int iTau, int iTauC, double dRhoFixed, double poisson_param, long initial_random_seed){
    FILE *fp;
    int i;

    Synapse * syn = (Synapse *)obj;

    fprintf(logfile, "Opening file for saving: %s\n", filename);
    fp = fopen(filename, "a");
    if (fp == NULL){
        perror("Error in saveSynapseOutputFile()");
    }
    else{
        fprintf(logfile, "Saving to file: %s\n", filename);

        fprintf(fp, "\n\n# Params:\n#    Cpre: %f, Cpost: %f, thetaD: %f, thetaP: %f, gammaD: %f, gammaP: %f, sigma: %f, Delay: %d, tau: %d, tauC: %d, rhoF: %f, poisson param: %f, seed: %ld\n", dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, iTau, iTauC, dRhoFixed, poisson_param, initial_random_seed);
        fprintf(fp,"#\n# Synaptic output, Synapse(%d):\n# t rho c preT postT\n", (*syn).ID);
        for (i = 0; i <= duration; i++){
            fprintf(fp, "%d %f %0.10lf %u %u\n", i, (*syn).rho[i], (*syn).c[i], (*syn).preT[i], (*syn).postT[i]);
        }
        fclose(fp);
        fprintf(logfile, "Completed saving\n");
    }

    return 0;
}


void loadSimulationParameters(int argc, char *argv[]){
    FILE* paramsfile;
    char paramName[FILE_NAME_LENGTH];
    double paramValue;

    strcpy(logfilearray, "output/logfile.log");
    logfilename = logfilearray;
    strcpy(outfilearray, "output/01_syn_%.3d.dat");
    outfilepattern = outfilearray;

    train_fn = train1;
    siT = 0;
    siID = 0;
    time_of_last_save = -1;

    if (argc < 2){
        printf("No command-line arguments passed....loading defaults.\n");
        printf("arg0 is %s\n", argv[0]);

        simulation_duration = 100;
        no_synapses = 1;

        initial_c = 0;
        initial_rho = 0;

        initial_random_seed = (-13);
        random_seed = initial_random_seed;

        iTau = 150000;  // measured in ms (equivalent to 150sec)
        iTauC = 20;  // measured in ms

        dRhoFixed = 0.5;

        dCpre = 1.0;
        dCpost = 4.0;
        dThetaD = 1.0;
        dThetaP = 2.0;
        dGammaD = 50;
        dGammaP = 100;

        dSigma = 3.50;
        iPreSpikeDelay = 19;
        poisson_param = 1.0/1000;
    }
    else{
        printf("arg[1] is %s....attempting to open it as parameters file.\n", argv[1]);
        paramsfile = fopen(argv[1], "r");
        if (paramsfile == NULL){
            perror("Error in loadSimulationParameters()");
        }
        else{
            printf("Loading from file: %s\n", argv[1]);
            while (!feof(paramsfile)) {
                if (fscanf(paramsfile, "%s %lf\n", paramName, &paramValue) != 2)
                    break;
                printf("Read: %s %f\n", paramName, paramValue);
                if (!strcmp(paramName, "SIMULATION_DURATION")){
                    simulation_duration = paramValue;
                }
                else if (!strcmp(paramName, "NO_SYNAPSES")){
                    no_synapses = paramValue;
                }
                else if (!strcmp(paramName, "INITIAL_C")){
                    initial_c = paramValue;
                }
                else if (!strcmp(paramName, "INITIAL_RHO")){
                    initial_rho = paramValue;
                }
                else if (!strcmp(paramName, "INITIAL_RANDOM_SEED")){
                    initial_random_seed = paramValue;
                    random_seed = paramValue;
                }
                else if (!strcmp(paramName, "TAU")){
                    iTau = (int) paramValue;
                }
                else if (!strcmp(paramName, "TAU_C")){
                    iTauC = (int) paramValue;
                }
                else if (!strcmp(paramName, "RHO_FIXED")){
                    dRhoFixed = paramValue;
                }
                else if (!strcmp(paramName, "C_PRE")){
                    dCpre = paramValue;
                }
                else if (!strcmp(paramName, "C_POST")){
                    dCpost = paramValue;
                }
                else if (!strcmp(paramName, "THETA_D")){
                    dThetaD = paramValue;
                }
                else if (!strcmp(paramName, "THETA_P")){
                    dThetaP = paramValue;
                }
                else if (!strcmp(paramName, "GAMMA_D")){
                    dGammaD = paramValue;
                }
                else if (!strcmp(paramName, "GAMMA_P")){
                    dGammaP = paramValue;
                }
                else if (!strcmp(paramName, "SIGMA")){
                    dSigma = paramValue;
                }
                else if (!strcmp(paramName, "PRE_SPIKE_DELAY")){
                    iPreSpikeDelay = (int) paramValue;
                }
                else if (!strcmp(paramName, "TRAIN_FUNCTION")){
                    if ((int) paramValue == 1){
                        train_fn = train1;
                    }
                    else if ((int) paramValue == 2){
                        train_fn = train2;
                    }
                    else if ((int) paramValue == 3){
                        train_fn = train3;
                    }
                    else if ((int) paramValue == 4){
                        train_fn = train4;
                    }
                    else if ((int) paramValue == 5){
                        train_fn = train5;
                    }

                    // - Pasted from (my) Pfister code

                    else if ((int) paramValue == 6){
                        train_fn = train6;
                    }
                    else if ((int) paramValue == 7){
                        train_fn = train7;
                    }
                    else if ((int) paramValue == 8){
                        train_fn = train8;
                    }
                    else if ((int) paramValue == 9){
                        train_fn = train9;
                    }
                    else if ((int) paramValue == 10){
                        train_fn = train10;
                    }
                    else if ((int) paramValue == 11){
                        train_fn = train11;
                    }
                    else if ((int) paramValue == 12){
                        train_fn = train12;
                    }
                    else if ((int) paramValue == 13){
                        train_fn = train13;
                    }
                    else if ((int) paramValue == 14){
                        train_fn = train14;
                    }
                    else if ((int) paramValue == 15){
                        train_fn = train15;
                    }
                    else if ((int) paramValue == 16){
                        train_fn = train16;
                    }
                    else if ((int) paramValue == 17){
                        train_fn = train17;
                    }
                    else if ((int) paramValue == 18){
                        train_fn = train18;
                    }
                    else if ((int) paramValue == 19){
                        train_fn = train19;
                    }
                    else if ((int) paramValue == 20){
                        train_fn = train20;
                    }
                    else if ((int) paramValue == 21){
                        train_fn = train21;
                    }
                    else if ((int) paramValue == 22){
                        train_fn = train22;
                    }
                    else if ((int) paramValue == 23){
                        train_fn = train23;
                    }
                    else if ((int) paramValue == 24){
                        train_fn = train24;
                    }
                    else if ((int) paramValue == 25){
                        train_fn = train25;
                    }

                    // - End of paste from (my) Pfister code

                    else if ((int) paramValue == 26){
                        train_fn = train26;
                    }
                    else if ((int) paramValue == 27){
                        train_fn = train27;
                    }

                    printf("DEBUG: train function loaded from param file\n");
                }
                else if (!strcmp(paramName, "POISSON_PARAM")){
                    poisson_param = paramValue;
                }
            }
            fclose(paramsfile);
        }

        if (argc > 2){
            printf("argv[2] is %s....attempting to use it as log file.\n", argv[2]);
            if (!strcmp(argv[2], "stdout")){
                strcpy(logfilearray, "stdout");
            }
            else{
                strcpy(logfilearray, "output/");
                strcat(logfilearray, argv[2]);
            }
            printf("Log file name: %s\n", logfilename);
        }
        if (argc > 3){
            printf("argv[3] is %s....attempting to use it as output file pattern.\n", argv[3]);
            strcpy(outfilearray, "output/");
            strcat(outfilearray, argv[3]);
            //outfilepattern = argv[3];
            printf("Output file pattern: %s\n", outfilepattern);
        }
    }

    // Make sure that directory 'output' exists
    if(mkdir("output",(S_IRUSR | S_IWUSR | S_IXUSR)) == -1){
        if (errno == EEXIST){
            printf("Directory 'output' already exists.\n");
        }
        else{
            perror("Error creating directory 'output'");
        }
    }
    // Open a log file for writing
    logfile = openLogFile(logfilename);
    fprintf(logfile, "\n\n---------------------\n  New Run\n---------------------\n");
}


int printToLog(FILE* fp, char* message){
    //printf("LOG: Saving message: %s\n", message);
    fprintf(fp, "%s", message);
    //printf("LOG: Completed saving\n");

    return 0;
}


FILE* openLogFile(char* filename){
    FILE *fp;
    if(!strcmp(filename, "stdout")){
        fp = stdout;
    }
    else{
        fp = fopen(filename, "a");
    }
    if (fp == NULL){
        perror("Error in openLogFile()");
    }
    else{
        printf("LOG: Log file opened for writing: %s\n", filename);
    }
    return fp;
}


int closeLogFile(FILE* fp){
    fclose(fp);
    return 0;
}

//// I'm only playing around here, was trying to
//// define object type as a parameter, turns out
//// I'd just be rewriting the basic libraries.
//int loadDataFile(char* filename, void *obj){
//    FILE *fp;
//
//    printf("Opening file for loading: %s\n", filename);
//    fp = fopen(filename, "r");
//    if (fp == NULL){
//        perror("Error in loadDataFile()");
//    }
//    else{
//        printf("Loading from file: %s\n", filename);
//        fscanf(fp, "%s\n", obj);
//        fclose(fp);
//    }
//
//    printf("Loaded data: %s", (char *)obj);
//
//    return 0;
//}
//
//// Note to self: when passing a (void *) argument, reassignment
//// to a local variable is unavoidable otherwise uncertain
//// compiler behaviour occurs (can't explain why)
//int saveOutputFile(char* filename, void *obj){
//    FILE *fp;
//    //double * local_obj = (double *)obj;
//    int * local_obj = (int *)obj;
//
//
//    printf("Opening file for saving: %s\n", filename);
//    fp = fopen(filename, "a");
//    if (fp == NULL){
//        perror("Error in saveOutputFile()");
//    }
//    else{
//        printf("Saving to file: %s\n", filename);
//        fprintf(fp,"Data value: %d\n", *local_obj);
//        fclose(fp);
//    }
//
//    // double
////    printf("output:%f\n", local_obj);
////    printf("value:%d\n", local_obj);
////    printf("local p value:%f\n", *local_obj);
////    printf("final value:%f\n", (double *)obj);
//
//    // int
//    printf("output:%d\n", local_obj);
//    printf("value:%d\n", local_obj);
//    printf("local p value:%d\n", *local_obj);
//    printf("final value:%d\n", (int *)obj);
//
//    return 0;
//}
