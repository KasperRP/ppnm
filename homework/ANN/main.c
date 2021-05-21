#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// Declare variables and functions from ann.c
typedef struct{
	int n;
	double (*f)(double);
	double(*fm)(double);
	double(*F)(double);
	double* params;
} ann;

ann* ann_alloc(int n, double (*f)(double), double (*fm)(double), double (*F)(double));
void ann_free(ann* network);
double ann_response(ann* network, double x);
double ann_der(ann* network, double x);
double ann_integ(ann* network, double x, double x0);
double cost_func(int d, double* p);
void ann_train(ann* network, int nx, double* xs, double* ys);

// Activation function and its derivative and its antiderivative (I choose x*exp(-x*x) among other possibilities)

double activation_func(double x){
	return x*exp(-x*x);
}

double activation_der(double x){
	return (1-2*x*x)*exp(-x*x);
}

double activation_integ(double x){
	return -exp(-x*x)*1./2;
}

// function to be created by interpolation

double func(double x){
	return sin(x)*exp(-x);
}
// derivative

double func_der(double x){
	return (cos(x)-sin(x))*exp(-x);
}
// antiderivative

double func_integ(double x, double x0){
	return -1./2*(sin(x)+cos(x))*exp(-x)+1./2*(sin(x0)+cos(x0))*exp(-x0);
}


int main(){
	// generate network
	int n=6; //number of neurons
	ann* network = ann_alloc(n,activation_func, activation_der, activation_integ);
	double xmin=0, xmax=2; // interval on x-axis
	// set parameters
	for(int i=0; i<network->n; i++){
		network->params[3*i]=xmin+(xmax-xmin)*i/(network->n-1);
		network->params[3*i+1]=1;
		network->params[3*i+2]=1;
	}

	// generate data points
	int nx = 50; // number of points
	double xs[nx];
	double ys[nx]; //tabulated points
	double yms[nx]; //derivative
	double Ys[nx]; //antiderivative
	for(int i=0; i<nx;i++){
		xs[i]=xmin+(xmax-xmin)*i/(nx-1);
		ys[i]=func(xs[i]);
		yms[i]=func_der(xs[i]);
		Ys[i]=func_integ(xs[i],xmin);
	}
	
	// train the network
	ann_train(network, nx, xs, ys);

	// Found optimized patameters
	FILE* Exc_AB = fopen("Exc_AB.txt", "w");
	for(int i=0; i<network->n; i++){
		double ai = network->params[3*i];
		double bi = network->params[3*i+1];
		double wi = network->params[3*i+2];
		fprintf(Exc_AB, "i=%i, ai, bi, wi = %g %g %g\n", i, ai, bi, wi);
	}

	// plot interpolation after training
	// points:
	FILE* points =fopen("points.txt", "w");
	for(int i=0; i<nx;i++) fprintf(points, "%g %g %g %g\n", xs[i], ys[i], yms[i], Ys[i]);

	// functions
	FILE* functions = fopen("functions.txt", "w");
	double dz=0.2;
	for(double z=xmin; z<=xmax; z+=dz){
		fprintf(functions, "%g %g %g %g\n", z, ann_response(network,z), ann_der(network,z), ann_integ(network,z,xmin));
	}
ann_free(network);
fclose(points);
fclose(functions);
fclose(Exc_AB);
return 0;
}

