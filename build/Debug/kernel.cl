/*
 *  kernel.cl
 *  XclNet
 *
 *  Created by David Higgins on 26/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


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


// Leaky integrate and fire kernel
//
__kernel void lif(
	__global float* input_v, // membrane voltage
	__global float* input_i, // input current
	__global float* input_gauss, // gaussian noise on membrane potential
	__global unsigned int* input_spike, // refractory period count down variable
	
	//__global float* output_v, // updated membrane voltage
	//__global unsigned int* output_spike, // refractory period countdown variable
	
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float r_m, // membrane resistance
	const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const float refrac_time, // duration of refractory period
	const float dt, // time step size
	const unsigned int no_lifs // number of lifs in simulation	
	)
{
	int i = get_global_id(0);
	if ( i < no_lifs ){
		float v = input_v[i];
		float input_current = input_i[i];
		float gauss = input_gauss[i];
		unsigned int time_since_spike = input_spike[i];
	
		float new_v;
		float dv = 0;
		float noise = 0;
		float tau_m = r_m * c_m;
	
//	if (time_since_spike == 0){
//		// A spike has just occurred, reset membrane voltage to reset potential
//		v =  v_reset;
//		// Apply leak current
//		dv = (-(v - v_rest) / (c_m * r_m));
//		if ( refrac_time == 0){
//			// No refractory period in operation so add input current
//			dv += (input_current / c_m);
//		}
//		// Apply noise
//		noise = sqrt(dt / tau_m) * sigma * gauss;
//		new_v = v + (dv * dt) + noise;
//		// Check lower threshold of membrane voltage
//		if (new_v < v_rest){
//			new_v = v_rest;
//		}		
//	}
//	else if (time_since_spike < refrac_time){
//		// Still in refractory period after last spike, but a spike has not just occurred
//		// v should have been reset on a previous time step
//		// Apply leak current
//		dv = (-(v - v_rest) / (c_m * r_m));
//		// Apply noise
//		noise = sqrt(dt / tau_m) * sigma * gauss;
//		new_v = v + (dv * dt) + noise;
//		// Apply lower threshold on membrane voltage
//		if (new_v < v_rest){
//			new_v = v_rest;
//		}
//	}
//	else{
//		// Not in refractory period, normal operation
//		// Apply leak current
//		dv = (-(v - v_rest) / (c_m * r_m));
//		// Apply input current
//		dv += (input_current / c_m);
//		// Apply noise
//		noise = sqrt(dt /tau_m) * sigma * gauss;
//		new_v = v + (dv * dt) + noise;
//		// Apply lower threshold to membrane voltage
//		if (new_v < v_rest){
//			new_v = v_rest;
//		}
//	}
	
		//TODO: decide initial value for time_since_spike, otherwise system always resets to V_reset upon initialisation
		if (time_since_spike == 0){
			// A spike has just occurred, reset membrane voltage to reset potential
			v = v_reset;
		}
	
		//CONSIDER: changed to >= to allow removal of logical OR which didn't work: (refrac_time==0)||
		if ( time_since_spike >= refrac_time ){
			// Apply leak current
			dv = (-(v - v_rest) / (c_m * r_m));
			// If refractory period is 0 OR if it's been longer than the refractory period since the last
			// spike, apply the external current
			dv += (input_current / c_m);
			// Apply noise
			//TODO: does the noise only apply during the refractory period?
			noise = sqrt(dt / tau_m) * sigma * gauss;
		}

		new_v = v + (dv * dt) + noise;
		// Apply lower threshold to membrane voltage
		if (new_v < v_rest){
			new_v = v_rest;
		}
	
		//Check if a spike has just occurred
		if (new_v > v_threshold){
			// A spike has just occurred, set time since last spike to 0
			time_since_spike = 0;
		}
		else{
			// No spike occurred, increment time since last spike
			time_since_spike++;
		}
	
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
	
	//TODO: change to RW buffers
	//__global float* output_rho, // updated synaptic efficacy
	//__global float* output_ca, // updated synaptic calcium concentration
	
	//TODO: opencl guarantees support for a minimum of 8 const args, need to work around this probably by passing a struct
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
		float gauss = input_gauss[i];
		unsigned int pre_spike = input_pre_spike[i];
		unsigned int post_spike = input_post_spike[i];
	
		float new_rho, new_ca;
		float drho = 0, dca = 0;
		//float noise = 0;
	
		unsigned int h_pot = 0;
		unsigned int h_dep = 0;
		unsigned int h_noise = 0;
	
		// Update Calcium concentration
		dca = (-ca/tau_ca) + (c_pre * pre_spike) + (c_post * post_spike);
		new_ca = ca + (dca * dt);
	
		// Update Synaptic efficacy
		if (new_ca > theta_p){
			h_pot = 1;
			h_noise = 1;
		}
		if (new_ca > theta_d){
			h_dep = 1;
			h_noise = 1;
		}
		//TODO: check if noise should also be divided by tau (I guess not, in which case I need to change this)
		drho = (-rho * (1.0 - rho) * (0.5 - rho)) + (gamma_p * (1.0 - rho) * h_pot) - (gamma_d * rho * h_dep) + (h_noise * sigma * sqrt(tau) * gauss);
		drho /= tau;
		new_rho = rho + (drho * dt);
	
		input_rho[i] = new_rho;
		input_ca[i] = new_ca;
	}
}
