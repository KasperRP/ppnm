#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

typedef struct{double x; int n;} besselpar;

double f(double t, void* params){
	besselpar par = *(besselpar*) params; // inspired by params.c in lecture 7-passf
	double x = (*par).x;
	int n = (*par).n;
	f=1/M_PI*cos(n*t-x*sin(t));
	return f;
}

double bessel(double x, int n) {
	gsl_function F;
	F.function = &f;
	F.params = &par;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc(limit);
	double a = 0, b = M_PI, acc = 1e-6, eps = 1e-6, result, error;
	gsl_integration_qag(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main() {
	double xmin = -2, xmax = 2;
	for(x=xmin;x<=xmax;<+=1.0/8);

