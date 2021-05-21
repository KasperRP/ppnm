#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
// Functions to solve linear system of equations
double dot(gsl_vector* x, gsl_vector* y);
double norm(gsl_vector* x);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);


void newton(
	void f(gsl_vector* x, gsl_vector* fx),
	gsl_vector* x,
	double eps){

int n = x -> size;
int LIMIT = 10000;
int steps = 0;

gsl_matrix* J = gsl_matrix_alloc(n,n);
gsl_matrix* R = gsl_matrix_alloc(n,n);
gsl_vector* Dx = gsl_vector_alloc(n); 
gsl_vector* fx = gsl_vector_alloc(n);
gsl_vector* df = gsl_vector_alloc(n);
gsl_vector* z = gsl_vector_alloc(n); 
gsl_vector* fz = gsl_vector_alloc(n);

while(steps<LIMIT){
f(x,fx); // For the given x, calculate f(x) and store in fx

// We create the Jacobian:.
for(int k=0; k<n; k++){
	double xk = gsl_vector_get(x,k);
	gsl_vector_set(x, k, xk+sqrt(DBL_EPSILON)); // Our finite difference is square root of machine epsilon
	f(x,df); // now df=f(x+sqrt())
        gsl_vector_sub(df,fx); // df=f(x+sqrt())-f(x), as it should	
	// Now we have what we need to calculate J
	for(int i=0; i<n; i++){
		double Jik = (gsl_vector_get(df,i)/sqrt(DBL_EPSILON)); // Finite difference
		gsl_matrix_set(J,i,k,Jik);}
	gsl_vector_set(x,k,xk);
}
// Next we solve the equation J*Dx=-f(x) using the routines from the homework "Linear Equations"
GS_decomp(J,R);
GS_solve(J, R, fx, Dx); // solution of J*Dx=f(x) now stored in Dx. We want the solution of J*Dx=-f(x):
gsl_vector_scale(Dx,-1.0);
// Then, we introduce lambda and use the conservative step lambda*Dx
double lambda = 1.0; // initial value
while(lambda>1./64){
gsl_vector_memcpy(z,x); // copy of x
gsl_vector_add(z,Dx); // z=x+Dx
f(z, fz); // f(x+Dx)
if( norm(fz)<(1-lambda/2)*norm(fx)) break;
lambda*=0.5;
gsl_vector_scale(Dx,0.5);
}

// Put final result into fx 
gsl_vector_memcpy(x,z);
gsl_vector_memcpy(fx,fz);
if(norm(Dx)<sqrt(DBL_EPSILON) || norm(fx)<eps ) break;
steps++;
}

gsl_matrix_free(J);
gsl_matrix_free(R);
gsl_vector_free(Dx);
gsl_vector_free(fx);
gsl_vector_free(df);
gsl_vector_free(z);
gsl_vector_free(fz);

}





