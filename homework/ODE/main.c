#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void rkstep23(
		void f(int n, double x, double* y, double* dydx), // f(x,y)=dy/dx. Solution stored in dydx
		int n, // order of equation / number of entries in vector y
		double x,
		double* yx, // y(x)
		double h, 
		double* yh, // y(x+h)
		double dy  // error
	     );

void driver(
		void f(int n, double x, double* y, double* dydx),
		int n,
		double a, // first point in interval
		double* ya, // y(a)
		double b, // end point
		double* yb, // y(b) - to be calculated
		double h, 
		double acc,
		double eps,
		char* path); // trajectory path 


// For harmonic oscillator

void f(int n, double x, double* y, double* dydx){
	dydx[0]=y[1];  // u''=-u corresponds to y1=y0' and y1=-y0' if y0=u and y1=u'
	dydx[1]=-y[0];
}

// For the SIR model
double N;
double T_c;
double T_r;

void fSIR(int n, double x, double* y, double* dydx){
	dydx[0]=-y[0]*y[1]/(N*T_c);
	dydx[1]=y[0]*y[1]/(N*T_c)-y[1]/T_r;
	dydx[2]=y[1]/T_r;
}

int main(){

	{	
	// For a harmonic oscillator solution over the interval [0,2*pi]:
	int n = 2;
	double a = 0;
	double b = 2*M_PI;
	double ya[n];
	double yb[n];
	double h = 0.001;
	double acc = 1e-4;
	double eps = 1e-4;
	char* path = "harmonic.txt";

	// initial values
	ya[0]=0;
	ya[1]=-1;

	driver(&f, n, a, ya, b, yb, h, acc, eps, path);
	}


	{

	// For the SIR model. I will plot different solutions with variating T_c
	
	int n = 3;
	double a = 0;
	double b = 50;
	double ya[n];
	double yb[n];
	double h = 0.001;
	double acc = 1e-4;
	double eps = 1e-4;
	char* path = "SIR1.txt";

	// initial values. Initially S should be comparable to the total population, I small, and R zero
	N=6e6; // Population of DK
	T_c=0.5;
	T_r=10;

	ya[0]=N;
	ya[1]=50;
	ya[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);


	path= "SIR2.txt";
	N=6e6;
	T_c=1;
	T_r=10;

	ya[0]=N;
	ya[1]=50;
	ya[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);

	path= "SIR3.txt";
	N=6e6;
	T_c=1.5;
	T_r=10;

	ya[0]=N;
	ya[1]=50;
	ya[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);

	FILE* comment_file = fopen("Comments_SIR.txt", "w");
	fprintf(comment_file, "An increasing contact time $T_c$ gives rise to the solutions being\nshifted along the time axis.\nAlso, the maximum infectious number decreases for increasing $T_c$, which makes sense.");	

	fclose(comment_file);
	}





	return 0;
}
