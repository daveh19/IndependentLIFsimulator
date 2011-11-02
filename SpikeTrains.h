// train1 is one pre- and no post- synaptic spikes
int train1(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// train2 is 20 pre-synaptic spikes and no post-synaptic spikes
int train2(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// train3 is Poisson(n) distributed pre-synaptic spiking
int train3(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// train4 is Poisson(n) distributed post-synaptic spiking
int train4(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// train5 is Poisson(n) distributed pre- and post-synaptic spiking
int train5(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// Dual spike shot noise simulation - train26 below


// - Pasted from (my) Pfister code

// Sjoestrom 2001
// f=0.1, dt=+10
int train6(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=0.1 dt=-10
int train7(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);

// f=10, dt=+10
int train8(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=10 dt=-10
int train9(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=20, dt=+10
int train10(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=20 dt=-10
int train11(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=40, dt=+10
int train12(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=40 dt=-10
int train13(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=50, dt=+10
int train14(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// f=50 dt=-10
int train15(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);

// Wang 2005 Triplets - Pre-Post-Pre
// t1=5, t2=5
int train16(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=10, t2=10
int train17(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=5, t2=15
int train18(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=15, t2=5
int train19(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// Wang 2005 Triplets - Post-Pre-Post
// t1=5, t2=5
int train20(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=10, t2=10
int train21(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=5, t2=15
int train22(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// t1=15, t2=5
int train23(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);

// Wang 2005 Quadruplets - Post-Pre-Pre-Post
// Tmin = 7, Tincrement=5
int train24(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// Wang 2005 Quadruplets - Pre-Post-Post-Pre
// Tmin = 7, Tincrement=5
int train25(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);

// - End pasted from Pfister code


// Dual spike shot noise simulation
// Pre-spike occurs as poisson process, post-spike occurs T timesteps after pre-spike
int train26(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
// Two shot noise processes, one with a single spike and one with a split 'dual' shape
// Single: pre-spike occurs as a poisson process, rate n1
// Dual: Pre-spike occurs as poisson process, post-spike occurs T timesteps after pre-spike
int train27(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration);
