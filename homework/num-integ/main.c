#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>


double quad(
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

// Functions to be integrated
int calls;
double g1(double x){
	calls++;
	return sqrt(x);}

double g2(double x){
	calls++;
	return 4*sqrt(1-x*x);}

double g3(double x){
	calls ++;
	return 1/sqrt(x);}

double g4(double x){
	calls ++;
	return log(x)/sqrt(x);}

// Writing functions again for GSL routine gsl_integration_qags

double g1_GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*sqrt(x);}

double g2_GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*4*sqrt(1-x*x);}

double g3_GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha/sqrt(x);}

double g4_GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*log(x)/sqrt(x);}



int main(){	 
	

	double a=0;
	double b=1;
	double acc=0.001;
	double eps=0.001;

	// For GSL routine
	int limit = 10000;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double alpha = 1.0;
	double result;
	double error;
	FILE* Exc_AB = fopen("Exc_AB.txt", "w");
	{
	calls=0;
	double Q=integrate(g1,a,b,acc,eps);
	fprintf(Exc_AB, "Integral of sqrt(x) from 0 to 1\n");
	fprintf(Exc_AB,"Analytical result = 2/3 \nMy numerical estimate = %g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	}
	
	{
	calls=0;
	gsl_function G;
	G.function = &g1_GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
	fprintf(Exc_AB,"GSL value = %f\nGSL number of calls = %d\n", result, calls);
	}

	{
	calls=0;
	double Q=integrate(g2,a,b,acc,eps);
	fprintf(Exc_AB,"\nIntegral of 4*sqrt(1-x^2) from 0 to 1\n");
	fprintf(Exc_AB,"Analytical result = pi = 3.1415926535897932384626433 \nMy numerical estimate = %0.25g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	}

	{
	calls=0;
	gsl_function G;
	G.function = &g2_GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
	fprintf(Exc_AB,"GSL value = %0.25f\nGSL number of calls = %d\n", result, calls);
	}	
	
	{
	calls=0;
	double Q=CCintegrate(g2,a,b,acc,eps);
	fprintf(Exc_AB,"Using Clenshaw-Curtis:\nMy numerical estimate = %0.25g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	fprintf(Exc_AB,"\nWe see that the Clenshaw-Curtis transform makes the result a little better than the adaptive routine, but not better than GSL\n");
	}
	
	{
	calls=0;
	double Q=integrate(g3,a,b,acc,eps);
	fprintf(Exc_AB,"\nIntegral of 1/sqrt(x) from 0 to 1\n");
	fprintf(Exc_AB,"Analytical result = 2 \nMy numerical estimate = %g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	}
	
	{
	calls=0;
	gsl_function G;
	G.function = &g3_GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
	fprintf(Exc_AB,"GSL value = %f\nGSL number of calls = %d\n", result, calls);
	}
	{
	calls=0;
	double Q=CCintegrate(g3,a,b,acc,eps);
	fprintf(Exc_AB,"Using Clenshaw-Curtis:\nMy numerical estimate = %g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	}

	{
	calls=0;
	double Q=integrate(g4,a,b,acc,eps);
	fprintf(Exc_AB,"\nIntegral of ln(x)/sqrt(x) from 0 to 1\n");
	fprintf(Exc_AB,"Analtical result = -4 \nMy numerical estimate = %g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps  = %g\n", acc, eps);
	}
	
	{
	calls=0;
	gsl_function G;
	G.function = &g4_GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
	fprintf(Exc_AB,"GSL value = %f\nGSL number of calls = %d\n", result, calls);
	}

	{
	calls=0;
	double Q=CCintegrate(g4,a,b,acc,eps);
	fprintf(Exc_AB,"Using Clenshaw-Curtis:\n");
	fprintf(Exc_AB,"My numerical estimate =%g\nNumber of calls = %d\n",Q,calls);
	fprintf(Exc_AB,"acc = %g\neps = %g\n", acc, eps);
	fprintf(Exc_AB,"\nFor the two latter integrals with singular integrands it is clear by the number of counts\nhow much more efficient the Clenshaw-Curtis transform makes the calculation.");
	}	

	gsl_integration_workspace_free(w);

	fclose(Exc_AB);
	return 0;
}
