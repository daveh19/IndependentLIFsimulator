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

