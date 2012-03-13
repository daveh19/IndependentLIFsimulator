#include <stdio.h>
#include <math.h>

#define PI (atan(1.)*4.)

#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define mwc ((znew<<16)+wnew)

#define shr3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define cong (jcong=69069*jcong+1234567)

#define kiss ((mwc^cong)+shr3)
#define uni (kiss*2.328306e-10)

typedef unsigned long UL;

static UL z=362436069, w=521288629, jsr=123456789, jcong=380116160;

typedef struct random_struct2{
	float value;
	unsigned int d_z;
	unsigned int d_w;
	unsigned int d_jsr;
	unsigned int d_jcong;
} Random2;
void my_two(Random2 *rdm){
	// Multiply-with-carry
	(*rdm).d_z = 36969*((*rdm).d_z&65535)+((*rdm).d_z>>16);
	(*rdm).d_w = 18000*((*rdm).d_w&65535)+((*rdm).d_w>>16);
	unsigned int d_mwc = ((*rdm).d_z<<16)+(*rdm).d_w;
	
	// 3-shift-register generator
	(*rdm).d_jsr^=((*rdm).d_jsr<<17);
	(*rdm).d_jsr^=((*rdm).d_jsr>>13);
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
void my_GetNormal(Random2 *rdm)
{
	// Use Box-Muller algorithm
	my_two(rdm);
	float u1 = (*rdm).value;
	my_two(rdm);
	float u2 = (*rdm).value;
	float r = sqrt( -2.0*log(u1) );
	float theta = 2.0*PI*u2;
	(*rdm).value = r*sin(theta);
}
/*
static unsigned long m_w = 521288629, m_z = 362436069;
unsigned int GetUint()
{
	m_z = 36969 * (m_z & 65535) + (m_z >> 16);
	m_w = 18000 * (m_w & 65535) + (m_w >> 16);
	return (m_z << 16) + m_w;
}
double GetUniform()
{
	// 0 <= u < 2^32
	unsigned int u = GetUint();
	// The magic number below is 1/(2^32 + 2).
	// The result is strictly between 0 and 1.
	return (u + 1.0) * 2.328306435454494e-10;
}
double GetNormal()
{
	// Use Box-Muller algorithm
	double u1 = GetUniform();
	double u2 = GetUniform();
	double r = sqrt( -2.0*log(u1) );
	double theta = 2.0*PI*u2;
	return r*sin(theta);
}*/


typedef struct random_struct{
	float value;
	unsigned int m_z;
	unsigned int m_w;
} RandomStruct;
// Multiply-with-Carry random number generator
// static unsigned int m_w = 521288629, m_z = 362436069;
unsigned int GetUint(RandomStruct *rnd)
{
	(*rnd).m_z = 36969 * ((*rnd).m_z & 65535) + ((*rnd).m_z >> 16);
	(*rnd).m_w = 18000 * ((*rnd).m_w & 65535) + ((*rnd).m_w >> 16);
	return ((*rnd).m_z << 16) + (*rnd).m_w;
}
// Uniform distribution in interval (0,1)
float GetUniform(RandomStruct *rnd)
{
	// 0 <= u < 2^32
	unsigned int u = GetUint(rnd);
	// The magic number below is 1/(2^32 + 2).
	// The result is strictly between 0 and 1.
	(*rnd).value = (u + 1.0) * 2.328306435454494e-10;
	return (*rnd).value;
}
// Gaussian(0,1) distribution using Box-Muller algorithm
float GetNormal(RandomStruct *rnd)
{
	// Use Box-Muller algorithm
	float u1 = GetUniform(rnd);
	float u2 = GetUniform(rnd);
	float r = sqrt( -2.0*log(u1) );
	float theta = 2.0*PI*u2;
	return r*sin(theta);
}

int main(void){
	RandomStruct rnd;
	rnd.m_w = 521288629;
	rnd.m_z = 362436069;
	Random2 rdm;
	rdm.d_z = 362436069;
	rdm.d_w = 521288629;
	rdm.d_jsr=123456789;
	rdm.d_jcong=380116160;
	for(int i = 0; i < 10; i++){
		my_two(&rdm);
		my_GetNormal(&rdm); // overwriting previous value
		printf("uni: %f, cook: %f, gaussian: %f, me_too: %f, my_gauss: %f\n", uni, GetUniform(&rnd), GetNormal(&rnd), rdm.value, rdm.value);
	}
	return 0;
}