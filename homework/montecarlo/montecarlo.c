#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#define RANDOM ((double)rand()/RAND_MAX)

void plainmc(
	int dim, // Number of dimensions	
	double* a, 
	double* b, 
	double f(int dim, double* x), 
	int N, // Number of points to be sampled
       	double* result, // pointer to variable where the result will be written
       	double* error){	// pointer to variable where the error will be written
	
	double V=1; //volume
	for(int i=0;i<dim;i++) V*=b[i]-a[i];
	
	double sum=0;
	double sum2=0;
	double x[dim];

	for(int i=0; i<N; i++){
	for(int i=0; i<dim;i++){
		x[i]=a[i]+RANDOM*(b[i]-a[i]);} // generates random x
		double fx=f(dim, x);
		sum+=fx;
		sum2+=fx*fx;}

	double mean = sum/N;
	double var = sum2/N-mean*mean;

	*result = mean*V;
	*error = sqrt(var/N)*V;
}

// Next I implement a Monte Carlo integrator that uses quasi-random sequences.
// I choose the Halton sequences (see notes):

// First implement van der corput sequence 

double corput(int n, int base){
	double q=0;
	double bk = (double) 1/base;
	while(n>0){
		q+=(n % base)*bk; // % means modulo
	       	n /= base;
	       	bk /= base;}
		return q;
	}
	
// Then two Halton sequences with different bases. We use the difference between these two
// as error estimation

double halton1(int n, int dim, int i){
	int base[]={2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67}; // The number must be coprimes. So we just take prime numbers.
	int maxdim = sizeof(base)/sizeof(int);
	assert(dim <= maxdim);
	return corput(n, base[i]);
}

double halton2(int n, int dim, int i){
	int base[]={71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163};
	int maxdim = sizeof(base)/sizeof(int);
	assert(dim <= maxdim);
	return corput(n, base[i]);
}

// Same Monte Carlo routine but now we use the above to generate points instead of the in-built random generator

void quasimc(
	int dim, 
	double* a,
	double* b,
	double f(int dim, double* x),
	int N,
	double* result,
	double* error){

	double V=1;
	for(int i=0; i<dim; i++) V*=b[i]-a[i];

	double sum=0;
	double x[dim];
	for(int i=0; i<N; i++){
		for(int j=0; j<dim; j++)  x[j]=a[j]+halton1(i, dim, j)*(b[j]-a[j]);

	 double fx = f(dim, x);
	 sum += fx;}

	double mean1 = sum/N;

	 sum=0;

	for(int i=0; i<N; i++){
		for(int j=0; j<dim; j++) x[j]=a[j]+halton2(i, dim, j)*(b[j]-a[j]);

	double fx = f(dim, x);
	sum += fx;}

	double mean2 = sum/N;

	*result = V*(mean1+mean2)/2; // We use the average of the two implementations
	*error = fabs(V*mean1-V*mean2); // Use absolute difference as error estimate
}
       		       

// Part C: recursive stratified sampling (see notes)

double strata(
	int dim, double f(int dim, double* x),
	double* a, double* b,
	double acc, double eps,
	int n_reuse, double mean_reuse)
{

int N = 16*dim;
double V = 1;
for(int i=0; i<dim; i++) V*=b[i]-a[i];
int nl[dim], nr[dim];
double x[dim];
double meanl[dim], meanr[dim];
double mean=0;

for(int i=0; i<dim; i++){
	meanl[i]=0;
	meanr[i]=0;
	nl[i]=0;
	nr[i]=0; }

for(int i=0; i<N; i++){
	for(int j=0; j<dim; j++) x[j]=a[j]+RANDOM*(b[j]-a[j]);
	double fx = f(dim,x);
	mean+=fx;

for(int j=0; j<dim; j++){
	if(x[j]>(a[j]+b[j])/2) { nr[j]++; meanr[j]+=fx; }
	else { nl[j]++; meanl[j]+=fx; }
}
}
mean /= N;
for(int i=0; i<dim; i++){
	meanl[i] /= nl[i];
	meanr[i] /= nr[i];
}

int idiv=0;
double maxvar=0;
for(int i=0; i<dim; i++){
	double var = fabs(meanr[i]-meanl[i]);
	if(var > maxvar){ maxvar=var; idiv=i; }
}

double integ = (mean*N+mean_reuse*n_reuse)/(N+n_reuse)*V;
double error = fabs(mean_reuse-mean)*V;
double tolerance = acc+eps*fabs(integ);
if(error < tolerance) return integ;

double a2[dim];
double b2[dim];
for(int i=0; i<dim; i++) a2[i]=a[i];
for(int i=0; i<dim; i++) b2[i]=b[i];
a2[idiv] = (a[idiv]+b[idiv])/2;
b2[idiv] = (a[idiv]+b[idiv])/2;

double integl = strata(dim, f, a, b2, acc/sqrt(2), eps, nl[idiv], meanl[idiv]);
double integr = strata(dim, f, a2, b, acc/sqrt(2), eps, nr[idiv], meanr[idiv]);
return integl+integr;
}



