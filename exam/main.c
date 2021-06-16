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


// simple function to be integrated
double Func1(double x, double y){
       return x*y; }

double d1(double x){
 	return 1;}

double u1(double x){
	return x;}

int main(){

double a = 0;
double b = 1;
double acc = 0.01;
double eps = 0.01;

FILE* testfile = fopen("testfile.txt", "w");
double I = integrate2D(Func1, a, b, d1, u1, acc, eps);
fprintf(testfile, "My estimate = %g\n", I);
fprintf(testfile, "Analytical result = -1/8");
fclose(testfile);
return 0;
}
