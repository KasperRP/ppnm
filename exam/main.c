#include<math.h>
#include<assert.h>
#include<stdio.h>

double integrate2D(
		double F(double, double),
		double a,
		double b,
		double d(double),
		double u(double),
		double acc,
		double eps);


// simple functions to be integrated
int calls;
double Func1(double x, double y){
	calls++;
       return x*y; } 

double Func2(double x, double y){
	calls++;
	return -x/(y*y); }

double Func3(double x, double y){
	calls++;
	return -1./(y*y); }

// Limit functions for the y-integral
double d1(double x){
 	return 1;}

double u1(double x){
	return x;}

double u2(double x){
	return sqrt(x);}

int main(){

double a = 0;
double b = 1;
double acc = 0.01;
double eps = 0.01;
FILE* testfile = fopen("testfile.txt", "w");
{
calls=0;
double I = integrate2D(Func1, a, b, d1, u1, acc, eps);
fprintf(testfile, "Integral of x*y where y runs from 1 to x and x from 0 to 1:\n");
fprintf(testfile, "My estimate = %g\n", I);
fprintf(testfile, "Analytical result = -1/8\n");
fprintf(testfile, "acc = %g\neps = %g\n", acc, eps);
fprintf(testfile, "Number of calls = %d\n", calls);
}

{
calls=0;
double I = integrate2D(Func2, a, b, d1, u1, acc, eps);
fprintf(testfile, "\nIntegral of -x/y^2 where y runs from 1 to x and x from 0 to 1:\n");
fprintf(testfile, "My estimate = %g\n", I);
fprintf(testfile, "Analytical result = 1/2\n");
fprintf(testfile, "acc = %g\neps = %g\n", acc, eps);
fprintf(testfile, "Number of calls = %d\n", calls);
}

{
calls=0;
double I = integrate2D(Func3, a, b, d1, u2, acc, eps);
fprintf(testfile, "\nIntegral of -1/y^2 where y runs from 1 to sqrt(x) and x from 0 to 1:\n");
fprintf(testfile, "My estimate = %g\n", I);
fprintf(testfile, "Analytical result = 1\n");
fprintf(testfile, "acc = %g\neps = %g\n", acc, eps);
fprintf(testfile, "Number of calls = %d\n", calls);
fprintf(testfile, "\nFor the last integral the number of calls is expected to be large\ncompared to the others, due to the presence of a singularity. For a faster\nconvergence, one could have included the Clenshaw-Curtis transformation in\nthe one-dimensional adaptive integrator.\n"); 
}

fclose(testfile);
return 0;
}
