#include "GeneralIncludes.h"
#include "ExternalCurrents.h"
#include "LIFNeuron.h"
#include "NumericalTools.h"


// current1 is a constant input current
int current1(double * Iext, unsigned int simulation_duration){
    int i;
    float mu = 0.5; //0.8; matches threshold
//    float sigma = 5.1; //CONSIDER: I think this should be 0.25 but that doesn't produce enough firing
//    float tau_m = dCm * dRm;
//    float gauss;

    for (i = 0; i < simulation_duration; i++){
        if( i > lif_time_of_last_save){
//            gauss = gasdev(&random_seed);
            Iext[i] = mu; // + (sigma * sqrt(tau_m) * gauss);
        }
    }

    return 0;
}

//// train2 is 20 pre-synaptic spikes and no post-synaptic spikes
//int train2(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
//    int i;
//    for (i = 0; i < 20; i++){
//        if( (i > lif_time_of_last_save) && (i < simulation_duration) ){
//            preT[i] = 1;
//            postT[i] = 0;
//        }
//    }
//    for (;i < simulation_duration; i++){
//        if( i > lif_time_of_last_save){
//            preT[i] = 0;
//            postT[i] = 0;
//        }
//    }
//
//    return 0;
//}
