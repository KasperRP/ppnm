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
		double b,	// end point
		//double* y,
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

double G;
double M1;
double M2;
double M3;
void fGRAV(int n, double x, double* y, double* dydx){
	// Vector y has 12 entries: position (x,y) and velocity (vx,vy) for the three different masses
	double x1=y[0]; double y1=y[1];
	double x2=y[2]; double y2=y[3]; // Not optimal notation
	double x3=y[4]; double y3=y[5];
	double vx1=y[6]; double vy1=y[7];
	double vx2=y[8]; double vy2=y[9];
	double vx3=y[10]; double vy3=y[11];

	double dx1dt=vx1; double dy1dt=vy1;
	double dx2dt=vx2; double dy2dt=vy2;
	double dx3dt=vx3; double dy3dt=vy3;

	double x12=x2-x1; double y12=y2-y1;
	double x13=x3-x1; double y13=y3-y1;
	double x23=x3-x2; double y23=y3-y2;

	double r12 = sqrt(x12*x12+y12*y12);
	double r13 = sqrt(x13*x13+y13*y13);
	double r23 = sqrt(x23*x23+y23*y23);

	double f12 = G*M1*M2/r12/r12;
	double f13 = G*M1*M3/r13/r13;
	double f23 = G*M2*M3/r23/r23;

	// Newtons Second Law:
	double dvx1dt = 1./M1*(f12*x12/r12+f13*x13/r13);
	double dvy1dt = 1./M1*(f12*y12/r12+f13*y13/r13);
	double dvx2dt = 1./M2*(-f12*x12/r12+f23*x23/r23);
	double dvy2dt = 1./M2*(-f12*y12/r12+f23*y23/r23);
	double dvx3dt = 1./M3*(-f13*x13/r13-f23*x23/r23);
	double dvy3dt = 1./M3*(-f13*y13/r13-f23*y23/r23);

	dydx[0]=dx1dt; dydx[1]=dy1dt;
	dydx[2]=dx2dt; dydx[3]=dy2dt;
	dydx[4]=dx3dt; dydx[5]=dy3dt;
	dydx[6]=dvx1dt; dydx[7]=dvy1dt;
	dydx[8]=dvx2dt; dydx[9]=dvy2dt;
	dydx[10]=dvx3dt; dydx[11]=dvy3dt;
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
	//double y[n];
	//y[0]=0; // initial values
	//y[1]=-1;

	driver(&f, n, a, ya, b, yb, h, acc, eps, path);
	//driver(f,n,a,b,y,h,acc,eps,path);
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
	//double y[n];
	//y[0]=N;
	//y[1]=50;
	//y[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);
	//driver(&fSIR, n, a, b, y, h, acc, eps, path);


	path= "SIR2.txt";
	N=6e6;
	T_c=1;
	T_r=10;

	ya[0]=N;
	ya[1]=50;
	ya[2]=0;
	//y[0]=N;
	//y[1]=50;
	//y[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);
	//driver(&fSIR, n, a, b, y, h, acc, eps, path);

	path= "SIR3.txt";
	N=6e6;
	T_c=1.5;
	T_r=10;

	ya[0]=N;
	ya[1]=50;
	ya[2]=0;
	//y[0]=N;
	//y[1]=50;
	//y[2]=0;

	driver(&fSIR, n, a, ya, b, yb, h, acc, eps, path);
	//driver(&fSIR, n, a, b, y, h, acc, eps, path);

	FILE* comment_file = fopen("Comments_SIR.txt", "w");
	fprintf(comment_file, "An increasing contact time $T_c$ gives rise to the solutions being\nshifted along the time axis.\nAlso, the maximum infectious number decreases for increasing $T_c$, which makes sense.");	

	fclose(comment_file);
	}

	{  // Gravitational Three-body problem
	int n = 12;
	double a=0;
	double ya[n];
	double b=6;
	double yb[n];
	double h = 0.1;
	double acc = 1e-3;
	double eps = 1e-3;	
	char* path = "GRAV.txt";

	G=1; M1=1; M2=1; M3=1;
	//double y[n];
	// Initial values
	ya[0]=-0.97000436; ya[1]=0.24308753; // x1 y1
	ya[2]=0; ya[3]=0; // x2 y2
	ya[4]=0.9700436; ya[5]=-0.24308753; // x3 y3
	ya[6]=0.4662036850; ya[7]=0.4323657300; // vx1 vy1
	ya[8]=-0.93240737; ya[9]=-0.86473146; // vx2 vy2
	ya[10]=0.4662036850; ya[11]=0.4323657300; // vx3 vy3

	driver(&fGRAV, n, a, ya, b, yb, h, acc, eps, path);
	//driver(&fGRAV, n, a, b, y, h, acc, eps, path);


	}



	return 0;
}
