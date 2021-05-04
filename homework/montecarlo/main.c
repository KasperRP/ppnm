#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void plainmc( int dim, double* a, double* b, double f(int dim, double* x), int N, double* result, double* error);

void quasimc( int dim, double* a, double* b, double f(int dim, double* x), int N, double* result, double* error);

// Interesting functions

double func1(int dim, double* r){
	assert(dim == 3);
	double x=r[0];
	double y=r[1];
	double z=r[2];
	return exp(x*y)*sin(M_PI*z);
}

double func2(int dim, double* r){
	assert(dim == 2 );
	double x=r[0];
	double y=r[1];
	return cos(x)*cos(x)*sin(y)*sin(y);
}

double func3(int dim, double* r){ // Stated in the exercise
	assert(dim == 3);
	double x=r[0];
	double y=r[1];
	double z=r[2];
	return 1/(1-cos(x)*cos(y)*cos(z))/(M_PI*M_PI*M_PI);
}

int main(){
	
	FILE* Exc_A = fopen("Exc_A.txt", "w");
	fprintf(Exc_A, "Testing the plain Monte Carlo integrator routine on different functions:\n");
	// Integral of func1
	{ 
	int dim = 3;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=1;
	}

	double result;
	double error;
	int N = 1e6;

	plainmc(dim, a, b, func1, N, &result, &error);
	fprintf(Exc_A, "\nIntegral of exp(x*y)*sin(pi*z) from [0,0,0] to [1,1,1]  = %g\n", result);
	fprintf(Exc_A, "Analytical result = 0.839003\n");
	fprintf(Exc_A, "Error = %g\nNumber of sampled points = %i\n", error,N);
	}
	
	// Integral of func2
	{
	int dim = 2;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=2*M_PI;
	}
	double result;
	double error;
	int N = 1e6;

	plainmc(dim, a, b, func2, N, &result, &error);
	fprintf(Exc_A, "\nIntegral of cos(x)^2*sin(y)^2 from [0,0] to [2pi, 2pi] = %g\n", result);
	fprintf(Exc_A, "Analytical result = pi^2\n");
	fprintf(Exc_A,"Error = %g\nNumber of sampled points = %i\n", error, N);
	}


	// Integral of func3 (The one stated in the exercise)
	{
	int dim = 3;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=M_PI;
	}
	double result;
	double error;
	int N = 1e6;

	plainmc(dim, a, b, func3, N, &result, &error);
	fprintf(Exc_A, "\nIntegral of 1/(pi^3*cos(x)*cos(y)*cos(z)) form [0,0,0] to [pi,pi,pi] = %0.25g\n", result);
	fprintf(Exc_A, "Stated Value = 1.3932039296856768591842462603255\n");
	fprintf(Exc_A, "Error = %g\nNumber of sampled points = %i\n", error, N);
	}
	
	FILE* Exc_B = fopen("Exc_B.txt", "w");
	fprintf(Exc_B, "Testing the quasi-random Monte Carlo routine on same functions as in A:\n");
	
	// Integral of func1
	{
	int dim = 3;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=1;
	}
	double result;
	double error;
	int N = 1e6;

	quasimc(dim, a, b, func1, N, &result, &error);
	fprintf(Exc_B, "\nIntegral of exp(x*y)*sin(pi*z) from [0,0,0] to [1,1,1] = %g\n", result);
	fprintf(Exc_B, "Analytical result = 0.839003\n");
	fprintf(Exc_B, "Error = %g\nNumber of sampled points = %i\n", error, N);
	}

	// Integral of func2
	{
	int dim = 2;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=2*M_PI;
	}
	double result;
	double error;
	int N = 1e6;

	quasimc(dim, a, b, func2, N, &result, &error);
	fprintf(Exc_B, "\nIntegral of cos(x)^2*sin(y)^2 from [0,0] to [2pi, 2pi] = %g\n", result);
	fprintf(Exc_B, "Analytical result = pi^2\n");
	fprintf(Exc_B, "Error = %g\nNumber of sampled points = %i\n", error, N);
	}
	
	// Testing error scaling of plain and quasi MC using my function func2
	
	FILE* err = fopen("err.txt", "w");
	{	
	int dim = 2;
	double a[dim];
	double b[dim];
	for(int i=0; i<dim; i++){
		a[i]=0;
		b[i]=2*M_PI;
	}
	double result_plain;
	double result_quasi;
	double error_plain;
	double error_quasi;
	for(int N=1000; N<100000; N+=1000){
		plainmc(dim, a, b, func2, N, &result_plain, &error_plain);
		quasimc(dim, a, b, func2, N, &result_quasi, &error_quasi);
		fprintf(err, "%i %g %g\n", N, error_plain, error_quasi);
	}

	}
	fclose(Exc_A);
	fclose(Exc_B);
	fclose(err);
	return 0;
	
	}
