#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void rkstep23(
		void f(int n, double x, double* y, double* dydx),
		int n, 
		double x,
		double* yx,
		double h,
		double* yh,
		double* dy);

void driver(
		void f(int n, double x, double* y, double* dydx),
		int n,
		double a,
		double* ya,
		double b,
		double* yb,
		double h,
		double acc,
		double eps,
		char* path);

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps);

static int ncalls;

// A simple one-dimensional function (A)

void f_one(gsl_vector* x, gsl_vector* fx){
	ncalls ++;
	double x0 = gsl_vector_get(x,0);
	double fx0 = x0*x0*x0-4*x0;
	gsl_vector_set(fx, 0, fx0);
}

// Two dimensional example (A)

void f_two(gsl_vector* x, gsl_vector* fx){
	ncalls++;
	double x0 = gsl_vector_get(x,0);
	double x1 = gsl_vector_get(x,1);
	double fx0 = x0*x0;
	double fx1 = x1*x1*x1;
	gsl_vector_set(fx, 0, fx0);
	gsl_vector_set(fx, 1, fx1);
}

// Gradient of Rosenbrock's valley function (A)


void Rosenbrock_grad(gsl_vector* r, gsl_vector* Rr){
	ncalls++;
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	double Rx = 2*(x-1)+400*(x*x-y)*x;
	double Ry = 200*(y-x*x);
	gsl_vector_set(Rr, 0, Rx);
	gsl_vector_set(Rr,1,Ry);
}

// Schrödinger ODE to be solved (B)
double e; // energy
void schroedinger(int n, double x, double* y, double* dydx){
	dydx[0]=y[1];
	dydx[1]=-2*e*y[0]+1.0/x*y[0];


	//dydx[0]=y[1];
	//dydx[1]=-y[0];
}

// Function giving M(e)=F_e(rmax) (B)
double rmax;
char* path;

void hydrogen(gsl_vector* x, gsl_vector* M){
	ncalls++;
	e = gsl_vector_get(x,0); // e is the variable of M
	int n = 2; // it's a second order ODE
	double a = 1e-3; // we avoid dividing by zero
	double b = rmax;

	//double a=0;
	//double b=rmax;
	
	double ya[n];
	double yb[n];
	double h = 0.001;
	double acc = 1e-4;
	double eps = 1e-4;

	// initial values - comes from f=r-r^2 and dfdr = 1-2r for small r
	ya[0]=a-a*a;
	ya[1]=1-2*a;
	

	//ya[0]=0;
	//ya[1]=-1;

	driver(&schroedinger, n, a, ya, b, yb, h, acc, eps, path);
	gsl_vector_set(M,0,yb[0]); // we find M as the solution evaluated at the endpoint.
	// Note that we only need the first component since we are not interested in the 
	// derivative of M.
}


int main(){

	{ // Part A
	FILE* Exc_A = fopen("Exc_A.txt", "w");
	fprintf(Exc_A, "We estimate the extremum of the Rosenbrock valley function\n");
	double eps = 0.01;
	int n = 2;
	gsl_vector* r = gsl_vector_alloc(n);
	double x0 = -2;
	double y0 = 8; // initial guess for roots
	gsl_vector_set(r,0, x0);
       gsl_vector_set(r, 1, y0);	
	
       ncalls = 0;
       newton(Rosenbrock_grad, r, eps);
       fprintf(Exc_A, "We start searching at (x,y) = (%g, %g)\n", x0, y0);
       fprintf(Exc_A, "The function is called %i times\n", ncalls);
       fprintf(Exc_A, "Found extremum of Rosenbrock's valley function:\n");
       for(int i=0; i<n; i++){
	       fprintf(Exc_A, "%g\n", gsl_vector_get(r,i));
       }
       fclose(Exc_A);
       gsl_vector_free(r);
	}

	{ // Part B
	rmax = 8;

	//rmax=2*M_PI;
	
	path = "hydrogen.txt";

	FILE* Exc_B = fopen("Exc_B.txt", "w");
	
	// We now have the solution F to 
	
	fclose(Exc_B);
	}

	return 0;
      
	}








     














