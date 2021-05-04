#include<math.h>
#include<assert.h>
#include<stdio.h>

double adaptint(
		double f(double x),
		double a,
		double b,
		double acc, 
		double eps,
		double f2,
		double f3,
		int nrec);

double integrate(
		double f(double x),
		double a,
		double b,
		double acc,
		double eps);

double CCintegrate(
		double f(double x),
		double a,
		double b,
		double acc,
		double eps);

int calls=0;
double g1(double x){
	calls++;
	return sqrt(x);};

double g2(double x){
	calls++;
	return 4*sqrt(1-x*x);};

double g3(double x){
	calls ++;
	return 1/sqrt(x);};


int main(){	// uses gcc nested functions as in the notes
	

	double a=0;
	double b=1;
	double acc=0.0001;
	double eps=0.0001;


	
	double intg1=integrate(g1,a,b,acc,eps);
	printf("Integral of sqrt(x) from 0 to 1\n");
	printf("Analytical result = 2/3 \nMy numerical estimate = %g\nNumber of calls = %d\n",intg1,calls);

	

	double intg2=integrate(g2,a,b,acc,eps);
	printf("\nIntegral of 4*sqrt(1-x^2) from 0 to 1\n");
	printf("Analytical result = pi \nMy numerical estimate = %g\nNumber of calls = %d\n",intg2,calls);
	
	double intg3=integrate(g3,a,b,acc,eps);
	printf("\nIntegral of 1/sqrt(x) from 0 to 1\n");
	printf("Analytical result = 2 \nMy numerical estimate = %g\nNumber of calls = %d\n",intg3,calls);

	double intg4=CCintegrate(g3,a,b,acc,eps);
	printf("\nIntegral of 1/sqrt(x) from 0 to 1 using Clenshaw-Curtis\n");
	printf("My numerical estimate =%g\nNumber of calls = %d\n",intg4,calls);

	return 0;
}
