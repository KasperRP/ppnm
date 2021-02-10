#include <math.h>
#include<tgmath.h>
#include <complex.h>
#include <stdio.h>
//#define M_E 2.7182818284590542354 //Tallene er allerede i de inkluderede biblioteker
//#define M_PI 3.1415926535897932384
int main(void) {
	printf("-------------------------------------------------------\n");
	double x;
	x=gamma(5);
	double y;
	y=j1(0.5);
	complex z;
	//z=csqrt(-2);
	//z=pow(M_E,I*M_PI);
	//z=pow(M_E,I);
	z=pow(I,M_E);
	//z=pow(I,I);
	printf("gamma(5)=%f\n",x);
	printf("j1(0.5)=%f\n",y);
	printf("Complex Number=%g +I%g\n",creal(z),cimag(z));
	printf("-------------------------------------------------------------------\n");
	double x_double = 1./9;
	float x_float = 1.f/9;
	long double x_long_double = 1.L/9;
	printf("%.25lg\n %.25g\n %.25Lg\n",x_double,x_float,x_long_double);

	return 0;
}
	
	
