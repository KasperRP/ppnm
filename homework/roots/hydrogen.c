#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>

void driver(
		void f(int n, double x, double* y, double* dydx),
		int n, double a, double* ya, double b, double* yb,
		double h, double acc, double eps, char* path);

double e;
void schroedinger(int n, double x, double* y, double* dydx){
	dydx[0]=y[1];
	dydx[1]=-2*e*y[0]+1.0/x*y[0];
}

double rmax = 8;
char* path = "hydrogen.txt";

void solver(gsl_vector* x, gsl_vector* M){
	e = gsl_vector_get(x,0);
	int n = 2;
	double a = 1e-3;
	double b = rmax;

	double ya[n];
	double yb[n];
	double h = 0.001;
	double acc = 1e-4;
	double eps = 1e-4;

	ya[0]=a-a*a;
	ya[1]=1-2*a;

	driver(&schroedinger, n , a, ya, b, yb, h, acc, eps, path);
	gsl_vector_set(M,0,yb[0]);
}
	
