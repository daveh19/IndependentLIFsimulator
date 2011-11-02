// Constructed from code in Numerical Recipes in C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <stddef.h>

void nrerror(char error_text[]);

float ran0(long *idum);
float ran1(long *idum);
float ran2(long *idum);
float ran3(long *idum);
float ran4(long *idum);
float gasdev(long *idum);
float expdev(long *idum);
float gamdev(int ia, long *idum);
float poidev(float xm, long *idum);
float bnldev(float pp, int n, long *idum);
