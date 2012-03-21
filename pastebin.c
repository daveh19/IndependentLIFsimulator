
//	(*syn_p).gamma_p = 725.085;
//	(*syn_p).gamma_d = 331.909;
//	(*syn_p).theta_p = 1.3;
//	(*syn_p).theta_d = 1.0;
//	//(*syn_p).delay = 4; // measured in multiples of dt
//	(*syn_p).sigma = 0; //3.35; //TODO: switch noise back on
//	(*syn_p).tau = 346.3615;
//	(*syn_p).tau_ca = 0.0226936;
//	(*syn_p).c_pre = 0.5617539;
//	(*syn_p).c_post = 1.23964;
//	(*syn_p).dt = 0.001;
//	//(*syn_p).no_syns = NO_SYNS;

/*
 // Generate links within Synapses
 // CONSIDER: this code could be much faster by using a temporarily much larger memory space and doing
 // processing withing generation of LIF links loops, then copying or resizing to real memory requirements
 printf("...links from synapses to pre and post synaptic neurons\n");
 (*syn_p).pre_lif = calloc(total_synapses, sizeof(signed int));
 (*syn_p).post_lif = calloc(total_synapses, sizeof(signed int));
 network_seed = NETWORK_SEED;
 int k = 0;
 for(int i = 0; i < NO_LIFS; i++){
 for(int j = 0; j < NO_LIFS; j++){
 if(i != j){
 if ((ran2(&network_seed)) < p){
 (*syn_p).pre_lif[k] = i;
 (*syn_p).post_lif[k] = j;
 k++;
 //printf("(%d)->(%d)\n", i, j);
 }
 }
 }
 }
 printf("k final (total synapses): %d\n", k);
 */	


//printf("Number of outgoing connections from lif(%d): %d\n", i, (*lif_p).no_outgoing_synapses[i]);
/*(*lif_p).outgoing_synapse_index[i] = calloc((*lif_p).no_outgoing_synapses[i], sizeof(signed int));
 for(int j = 0; j < k; j++){
 (*lif_p).outgoing_synapse_index[i][j] = out_links[j];
 //printf("Index(%d)(%d): %d \n", i, j, links[j]);
 }*/