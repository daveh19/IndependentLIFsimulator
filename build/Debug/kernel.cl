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
	
	__global float* output_v, // updated membrane voltage
	__global unsigned int* output_spike, // refractory period countdown variable
	
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float r_m, // membrane resistance
	const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const unsigned int refrac_time, // duration of refractory period
	const float dt // time step size	
	)
{
	int i = get_global_id(0);
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
	// Apply leak current
	dv = (-(v - v_rest) / (c_m * r_m));
	if ( (refrac_time == 0) || (time_since_spike > refrac_time) ){
		// If refractory period is 0 OR if it's been longer than the refractory period since the last
		// spike, apply the external current
		dv += (input_current / c_m);
	}
	// Apply noise
	noise = sqrt(dt / tau_m) * sigma * gauss;
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
	
	output_spike[i] = time_since_spike;
	output_v[i] = new_v;
}

