#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
typedef struct{ // structure containing the network
	int n; // number of neurons
	double(*f)(double); // activation function
	double(*fm)(double); // derivative of activation function 
	double(*F)(double); // antiderivative of activation function
	double* params; // parameters to be tuned during training
} ann;


ann* ann_alloc(int n, double(*f)(double),
		double(*fm)(double), double(*F)(double)
		){ // function to allocate memory of network
	ann* network = malloc(sizeof(ann)); // initialize with zeros
	network -> n = n;
	network -> f = f;
	network -> fm = fm;
	network -> F=F;
	network -> params = malloc(3*n*sizeof(double)); // three parameters (a,b,w) for each neuron
	return network;
}

void ann_free(ann* network){ // free allocated memory
	free(network -> params);
	free(network);
}


// Response function (Fp(x) in notes)
double ann_response(ann* network, double x){
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = network->params[3*i];
		double b = network->params[3*i+1];
		double w = network->params[3*i+2];
		sum+= network->f((x-a)/b)*w;
	}
	return sum;
}

// Derivative of responsfunction - note that we do not need this for tuning since we use the same parameters 
// for both the function, its derivative, and its antiderivative. But we still it defined to approximate derivative 
// of tabulated function.


double ann_der(ann* network, double x){
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = network->params[3*i];
		double b = network->params[3*i+1];
		double w = network->params[3*i+2];
		sum+= network->fm((x-a)/b)*w/b; // use some chain rule
	}
	return sum;
}

double ann_integ(ann* network, double x, double x0){ // antiderivative
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = network->params[3*i];
		double b = network->params[3*i+1];
		double w = network->params[3*i+2];
		sum+=network->F((x-a)/b)*w*b-network->F((x0-a)/b)*w*b;
	}
	return sum;
}



int qnewton(int d, double f(int, double*), double* x, double acc); // minimization function

// Cost function to be minimized in training function
static int N;
static double* X;
static double* Y;
static ann* NETWORK;
double cost_func(int d, double* p){ // avoiding nested functions
	assert(d==3*NETWORK->n);
	for(int i=0; i<d; i++)NETWORK->params[i]=p[i];
	double sum = 0;
	for(int k=0; k<N; k++){
		double fk = ann_response(NETWORK, X[k]);
		sum += (fk-Y[k])*(fk-Y[k]);
	}
	return sum/N;
}

void ann_train(ann* network, int nx, double* xs, double* ys){
	N=nx; X=xs; Y=ys;
	NETWORK=network;
	double acc = 1e-3;
	int d = 3*network->n;
	double p[d];
	for(int i=0;i<d;i++) p[i]=network->params[i];
	qnewton(d, cost_func, p, acc);
	for(int i=0;i<d;i++) network->params[i]=p[i];
}


