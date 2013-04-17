/*
 *  kernel.cl
 *  XclNet
 *
 *  Created by David Higgins on 26/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
 
#define PI (4.*atan(1.))

//#include "Random123/philox.h"

#include "Random123/philox.h"

/*typedef struct SynapseConstsStruct{
	float gamma_p;
	float gamma_d;
	float theta_p;
	float theta_d;
	unsigned int delay;
	float sigma;
	float tau;
	float tau_ca;
	float c_pre;
	float c_post;
	float dt;
	unsigned int no_syns;
} SynapseConsts;*/


// Simple compute kernel which computes the square of an input array
//
__kernel void square(                                  
   __global float* input,
   __global float* output,
   const unsigned int count)
{
	int i = get_global_id(0);
	if(i < count)
		output[i] = input[i] * input[i];
}


// Marsaglia's generators
// with corrected equations!
// taken from http://www.math.niu.edu/~rusin/known-math/99/RNG
typedef struct random_struct_marsaglia{
	float value;
	unsigned int d_z;
	unsigned int d_w;
	unsigned int d_jsr;
	unsigned int d_jcong;
} MarsagliaStruct;
// Get Uniform(0,1) pseudo-random value using technique published by Marsaglia
void Marsaglia_Uniform(MarsagliaStruct *rdm){
	// Multiply-with-carry
	(*rdm).d_z = 36969*((*rdm).d_z&65535)+((*rdm).d_z>>16);
	(*rdm).d_w = 18000*((*rdm).d_w&65535)+((*rdm).d_w>>16);
	unsigned int d_mwc = ((*rdm).d_z<<16)+(*rdm).d_w;
	
	// 3-shift-register generator
	(*rdm).d_jsr^=((*rdm).d_jsr<<13); //Fixed wrongly ordered shifts (see Kiss11)
	(*rdm).d_jsr^=((*rdm).d_jsr>>17);
	(*rdm).d_jsr^=((*rdm).d_jsr<<5);
	//unsigned int d_shr3 = (*rdm).d_jsr;
	
	// Congruential generator
	(*rdm).d_jcong = 69069*(*rdm).d_jcong+1234567;
	//unsigned int d_cong = (*rdm).d_jcong;
	
	// KISS generator (combines above three generators)
	unsigned int d_kiss = ((d_mwc^(*rdm).d_jcong)+(*rdm).d_jsr);
	
	// Convert to Uniform(0,1) distribution
	float d_uni = (d_kiss*2.328306e-10);
	
	(*rdm).value = d_uni;
}
// Convert Uniform(0,1) to Gaussian(0,1) using Box-Muller algorithm
void Marsaglia_GetNormal(MarsagliaStruct *rdm)
{
	// Use Box-Muller algorithm
	Marsaglia_Uniform(rdm);
	float r = sqrt( (float)(-2.0*log((*rdm).value)) );
	
	Marsaglia_Uniform(rdm); // Don't forget this second call!
	float theta = 2.0*PI*(*rdm).value;
	
	(*rdm).value = r*sin(theta);
}


typedef struct random_struct_mwc{
	float value;
	unsigned int m_z;
	unsigned int m_w;
} MWCRandomStruct;
// Multiply-with-Carry random number generator
// static unsigned int m_w = 521288629, m_z = 362436069;
// Code basically from http://www.codeproject.com/Articles/25172/Simple-Random-Number-Generation
unsigned int MWC_GetUint(MWCRandomStruct *rnd)
{
	(*rnd).m_z = 36969 * ((*rnd).m_z & 65535) + ((*rnd).m_z >> 16);
	(*rnd).m_w = 18000 * ((*rnd).m_w & 65535) + ((*rnd).m_w >> 16);
	return ((*rnd).m_z << 16) + (*rnd).m_w;
}
// Uniform distribution in interval (0,1)
float MWC_GetUniform(MWCRandomStruct *rnd)
{
	// 0 <= u < 2^32
	unsigned int u = MWC_GetUint(rnd);
	// The magic number below is 1/(2^32 + 2).
	// The result is strictly between 0 and 1.
	(*rnd).value = (u + 1.0) * 2.328306435454494e-10;
	return (*rnd).value;
}
// Gaussian(0,1) distribution using Box-Muller algorithm
void MWC_GetNormal(MWCRandomStruct *rnd)
{
	// Use Box-Muller algorithm
	float u1 = MWC_GetUniform(rnd);
	float u2 = MWC_GetUniform(rnd);
	float r = sqrt( (float)(-2.0*log(u1)) );
	float theta = 2.0*PI*u2;
	(*rnd).value = r*sin(theta);
}


// Leaky integrate and fire kernel
//
__kernel void lif(
	__global float* input_v, // membrane voltage
	__global float* input_i, // input current
	//__global float* input_gauss, // gaussian noise on membrane potential
	__global unsigned int* input_spike, // refractory period count up variable
	
	//State variables for random number generator
	/*__global unsigned int* d_z,
	__global unsigned int* d_w,
	__global unsigned int* d_jsr,
	__global unsigned int* d_jcong,*/
	
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float tau_m, // membrane time constant
	//const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const float refrac_time, // duration of refractory period
	const float dt, // time step size
	const unsigned int no_lifs, // number of lifs in simulation
	
	const unsigned int time_step, // used for indexing the random number generator
	const unsigned int random_seed, // seed for the random number generator
	
	//TODO: if gauss stream permanently removed then it should be removed from here, etc.
	__global float* output_gauss 
	)
	/*(
	__global float* input_v, // membrane voltage
	__global float* input_i, // input current
	__global float* input_gauss, // gaussian noise on membrane potential
	__global unsigned int* input_spike, // refractory period count down variable
	*/
	//State variables for random number generator
	/*__global unsigned int* d_z,
	__global unsigned int* d_w,
	__global unsigned int* d_jsr,
	__global unsigned int* d_jcong,*/
	/*
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float tau_m, // membrane resistance
	//const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const float refrac_time, // duration of refractory period
	const float dt, // time step size
	const unsigned int no_lifs, // number of lifs in simulation	
	
	const unsigned int time_step, // used for indexing the random number generator
	const unsigned int random_seed // seed for the random number generator
	)*/
{
	int i = get_global_id(0);
	if ( i < no_lifs ){
		float new_v;
		float dv;
		float noise;
		
		philox2x32_key_t key;
		philox2x32_ctr_t ctr;
		philox2x32_ctr_t rand_val;
		float2 uni_rand;
		
		float r, theta;
		float dave_temp;
		float my_random_value;
		
		float v = input_v[i];
		float input_current = input_i[i];
		unsigned int time_since_spike = input_spike[i];
		
		// Generate Gaussian(0,1) noise using Random123 library implementation		
		key.v[0] = i;
		ctr.v[0] = time_step;
		ctr.v[1] = random_seed;
		rand_val = philox2x32_R(10, ctr, key);
		// Convert to Uniform distribution (1/(2^32 +2))
		uni_rand.x = rand_val.v[0] * (1./(2^32 + 2));//2.328306435454494e-10;
		uni_rand.y = rand_val.v[1] * (1./(2^32 + 2));//2.328306435454494e-10;
		// Box-Muller transform
		r = sqrt( -2.0*log(uni_rand.x) );
		theta = 2.0 * PI * uni_rand.y;
		//my_random_value = r * sin(((float)theta);
		dave_temp = r * sin(theta);
		
		//dave_temp = dave_temp + 1;
		dv = 0;
		noise = 0;
		//my_random_value = 1.34;
		
		//my_random_value = my_random_value + 3.;
		
		// Generate Gaussian(0,1) noise
		/*MarsagliaStruct rnd;
		rnd.d_z = d_z[i];
		rnd.d_w = d_w[i];
		rnd.d_jsr = d_jsr[i];
		rnd.d_jcong = d_jcong[i];
		Marsaglia_GetNormal(&rnd);*/
	
		
		//float tau_m = r_m * c_m;
	
		//REMINDER: initialise time_since_spike to refrac_time in main program,
		// otherwise system always resets to V_reset upon initialisation
		if (time_since_spike == 0){
			// A spike has just occurred, reset membrane voltage to reset potential
			v = v_reset;
		}
	
		// If refractory period is 0 OR if it's been longer than the refractory period since the last spike
		//CONSIDER: changed to >= to allow removal of logical OR which didn't work: (refrac_time==0)||
		if ( time_since_spike >= refrac_time ){
			// Apply leak current
			dv = (-(v - v_rest) / tau_m);
			// Apply the external current
			// Note: I use one input current variable (to cut down on streams to GPU)
			//  an external current/voltage should be added directly to this variable (outside the kernel)
			//  a synaptic current/voltage step should be multiplied by (tau_m/dt), for a delta spike, before adding to this variable,
			//  in order to counter rescaling which happens on next three lines of executable code.
			// input_current is treated as a voltage step, despite the variable name, hence the division by tau_m
			dv += (input_current / tau_m);
			// Apply noise
			//noise = sqrt(dt / tau_m) * sigma * rnd.value;
			noise = sqrt(dt / tau_m) * sigma * my_random_value;
		}

		new_v = v + (dv * dt) + noise;
		// Apply lower threshold to membrane voltage (no longer desired)
		/*if (new_v < v_rest){
			new_v = v_rest;
		}*/
	
		//Check if a spike has just occurred
		if (new_v > v_threshold){
			// A spike has just occurred, set time since last spike to 0
			time_since_spike = 0;
		}
		else{
			// No spike occurred, increment time since last spike
			time_since_spike++;
		}
		
		/*d_z[i] = rnd.d_z;
		d_w[i] = rnd.d_w;
		d_jsr[i] = rnd.d_jsr;
		d_jcong[i] = rnd.d_jcong;*/
		
		output_gauss[i] = dave_temp;

		input_spike[i] = time_since_spike;
		input_v[i] = new_v;
	}
}


// Graupner 2012 Synapse kernel
//
__kernel void synapse( 
	__global float* input_rho, // synaptic efficacy
	__global float* input_ca, // synaptic calcium concentration
	__global float* input_gauss, // gaussian noise on synaptic efficacy
	__global unsigned int* input_pre_spike, // number of pre-synaptic spikes occurring at time t-D
	__global unsigned int* input_post_spike, // number of post-synaptic spikes at time t

	//State variables for random number generator
	__global unsigned int* d_z,
	__global unsigned int* d_w,
	__global unsigned int* d_jsr,
	__global unsigned int* d_jcong,
	
	//const SynapseConsts *syn_const, // struct of constants required for synapse
	//const unsigned int no_syns
	
	//TODO: opencl guarantees support for a minimum of only 8 const args, probably safer to pass params by passing a struct
	const float gamma_p, // potentiation learning rate
	const float gamma_d, // depression learning rate
	const float theta_p, // potentiation threshold
	const float theta_d, // depression threshold
	const float tau, // time constant for synaptic efficacy
	const float tau_ca, // time constant for calcium concentration
	const float c_pre, // increase in Ca by pre-synaptic spike
	const float c_post, // increase in Ca by post-synaptic spike
	const float sigma, // size of noise
	const float dt, // time step size	
	const unsigned int no_syns // number of synapses in simulation
	)
{
	int i = get_global_id(0);
	
	if (i < no_syns){
		float rho = input_rho[i];
		float ca = input_ca[i];
		unsigned int pre_spike = input_pre_spike[i];
		unsigned int post_spike = input_post_spike[i];
		
		// Generate Gaussian(0,1) noise
		MarsagliaStruct rnd;
		rnd.d_z = d_z[i];
		rnd.d_w = d_w[i];
		rnd.d_jsr = d_jsr[i];
		rnd.d_jcong = d_jcong[i];
		Marsaglia_GetNormal(&rnd);
	
		float new_rho, new_ca;
		float drho = 0., dca = 0.;
		float noise = 0.;
	
		unsigned int h_pot = 0;
		unsigned int h_dep = 0;
		unsigned int h_noise = 0;
	
		// Update Calcium concentration
		//TODO: surely the effects of c_pre and c_post should not be dependent on dt
		dca = (-ca/tau_ca); // + (c_pre * pre_spike) + (c_post * post_spike);
		new_ca = ca + (dca * dt) + (c_pre * pre_spike) + (c_post * post_spike);
	
		// Update Synaptic efficacy
		if (new_ca > theta_p){
			h_pot = 1;
			h_noise = 1;
		}
		if (new_ca > theta_d){
			h_dep = 1;
			h_noise = 1;
		}

		// Calculate rho update
		drho = (-rho * (1.0 - rho) * (0.5 - rho)) + (gamma_p * (1.0 - rho) * h_pot) - (gamma_d * rho * h_dep);
		drho /= tau;
		// Calculate noise
		noise = (h_noise * sigma * sqrt(dt/tau) * rnd.value);
		// Calculate new rho value
		new_rho = rho + (drho * dt) + noise;
	
		// Zeroing arrays preT and postT (for use in main)
		input_pre_spike[i] = 0;
		input_post_spike[i] = 0;
		
		// Final output
		d_z[i] = rnd.d_z;
		d_w[i] = rnd.d_w;
		d_jsr[i] = rnd.d_jsr;
		d_jcong[i] = rnd.d_jcong;
		input_gauss[i] = rnd.value;
		
		//input_rho[i] = new_rho;
		//TODO: re-enable synaptic weight change
		input_rho[i] = 1;
		input_ca[i] = new_ca;
		//TODO: double check numerical output for a given pre_spike or post_spike
	}
}
