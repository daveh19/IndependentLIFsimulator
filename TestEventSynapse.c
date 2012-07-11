/*
 *  TestEventSynapse.c
 *  XclNet
 *
 *  Created by David Higgins on 11/07/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralIncludes.h"
#include "cl_LIFNeuron.h"
#include "cl_Synapse.h"
//#include "HandleOpenCL.h"
#include "NumericalTools.h"

#include "DataReporters.h"

void updateEventBasedSynapse(cl_Synapse *syn, SynapseConsts *syn_const, int syn_id, int current_time);


int main(void){
	printf("Begin\n");
	
	long uniform_synaptic_seed = UNIFORM_SYNAPTIC_SEED;
	
	cl_Synapse syn;
	cl_Synapse *syn_p = &syn;
	SynapseConsts syn_const;
	SynapseConsts *syn_const_p = &syn_const;
	
	(*syn_const_p).no_syns = 1;
	
	(*syn_p).rho = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).rho_initial = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).ca = malloc(sizeof(float) * (*syn_const_p).no_syns);
	(*syn_p).gauss = calloc((*syn_const_p).no_syns, sizeof(float));
	(*syn_const_p).delay = SYN_CALCIUM_DELAY; // measured in multiples of dt
	
	(*syn_p).time_of_last_update = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	(*syn_p).preT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	(*syn_p).postT = calloc((*syn_const_p).no_syns, sizeof(unsigned int));
	
	(*syn_const_p).gamma_p = SYN_GAMMA_P;
	(*syn_const_p).gamma_d = SYN_GAMMA_D;
	(*syn_const_p).theta_p = SYN_THETA_P;
	(*syn_const_p).theta_d = SYN_THETA_D;
	(*syn_const_p).sigma = SYN_SIGMA;
	(*syn_const_p).tau = SYN_TAU;
	(*syn_const_p).tau_ca = SYN_TAU_CA;
	(*syn_const_p).c_pre = SYN_C_PRE;
	(*syn_const_p).c_post = SYN_C_POST;
	(*syn_const_p).dt = SYN_DT;
	
	for(int i = 0; i < (*syn_const_p).no_syns; i++){
		//(*syn_p).rho[i] = SYN_RHO_INITIAL;
		(*syn_p).rho[i] = (*syn_p).rho_initial[i] = 0.17470929; //ran2(&uniform_synaptic_seed);
		
		(*syn_p).ca[i] = 1.4; //SYN_CA_INITIAL;
		/*(*rnd_syn_p).d_z[i] = 362436069 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_w[i] = 521288629 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_jsr[i] = 123456789 - i + PARALLEL_SEED;
		(*rnd_syn_p).d_jcong[i] = 380116160 - i + PARALLEL_SEED;*/
	}
	
	(*syn_p).preT[0] = 1;
	(*syn_p).postT[0] = 0;
	int t = 0;
	printf("Before, t: %d, rho: %f, ca: %f\n", t, (*syn_p).rho[0], (*syn_p).ca[0]);
	t = 1644;
	updateEventBasedSynapse(syn_p, syn_const_p, 0, t);
	printf("After, t: %d, rho: %f, ca: %f\n", t, (*syn_p).rho[0], (*syn_p).ca[0]);
	
	printf("Done\n");
	
	return 0;
}



void updateEventBasedSynapse(cl_Synapse *syn, SynapseConsts *syn_const, int syn_id, int current_time){
	static long gaussian_synaptic_seed = GAUSSIAN_SYNAPTIC_SEED;
	float theta_upper = fmax((*syn_const).theta_d, (*syn_const).theta_p);
	float theta_lower = fmin((*syn_const).theta_d, (*syn_const).theta_p);
	float gamma_upper = fmax((*syn_const).gamma_d, (*syn_const).gamma_p);
	float gamma_lower = fmin((*syn_const).gamma_d, (*syn_const).gamma_p);
	/*float theta_upper = (*syn_const).theta_p;
	 float theta_lower = (*syn_const).theta_d;
	 float gamma_upper = (*syn_const).theta_p;
	 float gamma_lower = (*syn_const).theta_d;*/
	float w_stoch, w_deter, w;
	float c_initial, c_end;
	
	float time_since_update = (*syn_const).dt * (current_time - (*syn).time_of_last_update[syn_id]);
	
	c_initial = (*syn).ca[syn_id];
	w = (*syn).rho[syn_id];
	w_stoch = w_deter = 0;
	
	//if(syn_id == RECORDER_SYNAPSE_ID){
	printf("(SYN %d) seed: %ld, w_initial: %f, c_initial: %f, ", syn_id, gaussian_synaptic_seed, (*syn).rho[syn_id], c_initial);
	/*if(time_since_update > (*syn_const).dt){ // for graphing, fill in Ca value just before potential Ca influx
		c_end = c_initial * exp(-((double)(time_since_update - (*syn_const).dt) / (*syn_const).tau_ca));
		//TODO: print this out its the Recorder Synapse
		printf("time_since_update: %f, c_end before influx: %f, ", time_since_update, c_end);
		//(*syn).ca[current_time - 1] = c_end;
	}*/
	//}
	c_end = c_initial * exp(-((double)(time_since_update) / (*syn_const).tau_ca));
	printf("time_since_update: %f, c_end before influx: %f, ", time_since_update, c_end);
	
	//CONSIDER: test for time_since_update > 0 for rest of function (probably would take more clock cycles than allowing the calculation to proceed on that rare occasion)
	float t_upper, t_lower, t_deter;
	if (c_initial > theta_upper){
		if(c_end > theta_upper){
			//update tupper, tlower, tdeter and call stochastic update
			t_upper = time_since_update;
			t_lower = 0;
			t_deter = 0;
		}
		else if (c_end > theta_lower){ // && c_end <= theta_upper
			//update tupper, tlower, tdeter and call stochastic update
			t_upper = (*syn_const).tau_ca * log( c_initial/theta_upper );
			t_lower = time_since_update - t_upper;
			t_deter = 0;
		}
		else{ // c_end <= theta_lower
			//update tupper, tlower, tdeter and call stochastic update, then call deterministic update
			t_upper = (*syn_const).tau_ca * log( c_initial/theta_upper );
			t_lower = (*syn_const).tau_ca * log( theta_upper/theta_lower );
			t_deter = time_since_update - t_upper - t_lower;
		}
	}
	else if (c_initial <= theta_lower){
		//update tupper=0, tlower=0, tdeter and call deterministic update
		t_upper = 0;
		t_lower = 0;
		t_deter = time_since_update;
	}
	else if (c_end <= theta_lower){ // && c_initial > theta_lower && c_initial <= theta_upper
		//update tupper, tlower, tdeter and call stochastic update, then call deterministic update
		t_upper = 0;
		t_lower = (*syn_const).tau_ca * log( c_initial/theta_lower );
		t_deter = time_since_update - t_lower;
	}
	else{ // c_initial > theta_lower && c_initial <= theta_upper && c_end > theta_lower && c_end <= theta_upper
		//update tupper, tlower, tdeter and call stochastic update
		t_upper = 0;
		t_lower = time_since_update;
		t_deter = 0;
	}
	
	// Weight update
	/*if(t_lower > 0 || t_upper > 0){
	 float GammaP, GammaD, t_b, w_bar, tau_prime, sig_bar, sig_sq;
	 // Lower threshold depression, upper threshold potentiation
	 GammaP = (t_upper) * (*syn_const).gamma_p;
	 GammaD = (t_upper + t_lower) * (*syn_const).gamma_d;
	 t_b = t_upper + t_lower;
	 
	 w_bar = GammaP / (GammaD + GammaP);
	 tau_prime = (*syn_const).tau / (GammaD + GammaP);
	 w_mean = w_bar + (w - w_bar) * exp(-t_b/tau_prime);
	 
	 sig_bar = ((*syn_const).sigma / (2 * (GammaD + GammaP) ) );
	 sig_sq = pow(sig_bar,2) * (1 - exp(-(2*t_b)/tau_prime));
	 w_stoch = 0; // gaussian(0,sig_sq) distribution 
	 w = w_mean + w_stoch; // update here so deterministic update can follow on from stochastic one
	 }*/
	// Stochastic update
	double rnd;
	if (t_upper > 0){
		float rho_bar, in_exp, random_part, my_exp;
		rho_bar = (gamma_upper / (gamma_lower + gamma_upper));
		in_exp = -(t_upper * (gamma_lower + gamma_upper)) / (*syn_const).tau;
		my_exp = exp(in_exp);
		
		w_stoch = (gamma_upper / (gamma_lower + gamma_upper)) * ( 1 - my_exp);
		w_stoch += w * my_exp;
		printf("\nt_upper: %f, rho_bar: %f, in_exp: %f, my_exp: %f, w_stoch: %f, ", t_upper, rho_bar, in_exp, my_exp, w_stoch);
		rnd = gasdev(&gaussian_synaptic_seed);
		//rnd = gasdev(&gaussian_synaptic_seed);
		printf("rnd1: %f, ", rnd);
		random_part = (*syn_const).sigma * rnd * sqrt( (1 - exp(-(2 * (gamma_lower + gamma_upper) * t_upper) / (*syn_const).tau) ) / ( (2 * (gamma_lower + gamma_upper) ) ) );
		w_stoch += (*syn_const).sigma * rnd * sqrt( (1 - exp(-(2 * (gamma_lower + gamma_upper) * t_upper) / (*syn_const).tau) ) / ( (2 * (gamma_lower + gamma_upper) ) ) );
		printf("random: %f, w_stoch: %f\n", random_part, w_stoch);
		w = w_stoch;
	}
	if (t_lower > 0){
		float in_exp, my_exp, random_part;
		in_exp = -(t_lower * gamma_lower) / (*syn_const).tau;
		my_exp = exp(in_exp);
		printf("\nt_lower: %f, w: %f, in_exp: %f, my_exp: %f, ", t_lower, w, in_exp, my_exp);
		w_stoch = w * exp(-(t_lower * gamma_lower) / (*syn_const).tau);
		printf("w_stoch: %f, ", w_stoch);
		rnd = gasdev(&gaussian_synaptic_seed);
		printf("rnd2: %f, ", rnd);
		random_part = (*syn_const).sigma * rnd * sqrt( (1 - exp(-(2 * gamma_lower * t_lower) / (*syn_const).tau) ) / (2 * gamma_lower) );
		w_stoch += (*syn_const).sigma * rnd * sqrt( (1 - exp(-(2 * gamma_lower * t_lower) / (*syn_const).tau) ) / (2 * gamma_lower) );
		printf("random: %f, w_stoch: %f\n", random_part, w_stoch);
		w = w_stoch;
	}
	// Deterministic update
	if (t_deter > 0){
		float denominator = w * (w - 1);
		float numerator = pow(w - 0.5, 2);
		float X_0 = numerator / denominator;
		printf("\nt_deter: %f, w: %f, num: %f, den: %f, X_0: %f\n", t_deter, w, numerator, denominator, X_0);
		float in_exp, my_exp, X_exp, denom2, division, in_sqt, my_sqt, multiple, whole;
		in_exp = t_deter/(2 * (*syn_const).tau);
		my_exp = exp( in_exp );
		X_exp = X_0 * my_exp;
		denom2 = (X_exp - 1.);
		division = 1 / denom2;
		in_sqt = (1. + division ) ;
		my_sqt = sqrt(in_sqt);
		multiple = 0.5 * my_sqt;
		if (w < 0.5){
			whole = 0.5 - multiple;
			w_deter = 0.5 - (0.5 * sqrt( (1. + (1. / (X_0 * exp( t_deter/(2 * (*syn_const).tau) ) - 1.)) ) ) );
		}
		else{
			whole = 0.5 + multiple;
			w_deter = 0.5 + (0.5 * sqrt( (1. + (1. / (X_0 * exp( t_deter/(2 * (*syn_const).tau) ) - 1.)) ) ) );
		}
		printf("in_exp: %f, my_exp: %f, X_exp: %f, denom2: %f, division: %f, in_sqt: %f, my_sqt: %f, multiple: %f, w_deter: %f\n", in_exp, my_exp, X_exp, denom2, division, in_sqt, my_sqt, multiple, whole);
		w = w_deter;
	}
	
	c_end = c_end + ((*syn).preT[syn_id] * (*syn_const).c_pre) + ((*syn).postT[syn_id] * (*syn_const).c_post);
	//if(syn_id == RECORDER_SYNAPSE_ID){
	printf("after influx: %f, w_final: %f\n", c_end, w);
	//}
	// Reset preT and postT, so that calcium influx can only be applied once!
	(*syn).preT[syn_id] = 0;
	(*syn).postT[syn_id] = 0;
	(*syn).time_of_last_update[syn_id] = current_time;
	(*syn).ca[syn_id] = c_end;
	//TODO: should I put hard bounds on rho? (sigma=3.35 is too large otherwise)
	(*syn).rho[syn_id] = w;
}