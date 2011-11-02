#include "GeneralIncludes.h"
#include "NetworkTopology.h"
//#include "LIFNeuron.h"
#include "NumericalTools.h"


// topology 0 is a totally disconnected network (for simulator testing)
int topology0(signed int **adj, int no_neurons, long random_seed){
    int i, j;

    for (i = 0; i < no_neurons; i++){
        for (j = 0; j < no_neurons; j++){
            adj[i][j] = -1;
        }
    }

    return 0;
}


// topology 1 is an all-to-all connected network, including self-connections
int topology1(signed int **adj, int no_neurons, long random_seed){
    int i, j;
    int index = 0;

    for (i = 0; i < no_neurons; i++){
        for (j = 0; j < no_neurons; j++){
            //printf("%d ", adj[i][j]);
            adj[i][j] = index;
            index++;
        }
        //printf("\n");
    }

    return index;
}


// topology 2 is an all-to-all connected network, with no self-connections
int topology2(signed int **adj, int no_neurons, long random_seed){
    int i, j;
    int index = 0;

    for (i = 0; i < no_neurons; i++){
        for (j = 0; j < no_neurons; j++){
            //printf("%d ", adj[i][j]);
            if (i != j){
                adj[i][j] = index;
                index++;
            }
            else{
                adj[i][j] = -1;
            }
        }
        //printf("\n");
    }

    return index;
}


// topology 10 : the probability of a connection from a->b is 0.1, no self connections
int topology10(signed int **adj, int no_neurons, long random_seed){
    int i, j;
    int index = 0;
    float ran;

    for (i = 0; i < no_neurons; i++){
        for (j = 0; j < no_neurons; j++){
            if (i != j){
                ran = ran2(&random_seed);
                if ((ran) < 0.100000001){
                    //printf("DEBUG: ran (a) %f\n", ran);
                    adj[i][j] = index;
                    index++;
                }
                else{
                    //printf("DEBUG: ran (b) %f\n", ran);
                    adj[i][j] = -1;
                }
            }
            else{
                //printf("DEBUG: on diagonal, no self connections\n");
                adj[i][j] = -1;
            }
        }
    }

    return index;
}
