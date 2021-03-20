#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include"binsearch.h"

typedef struct {gsl_vector* x, *y, *b, *c, *d;} cinterp;

cinterp* cinterp_alloc(gsl_vector* x, gsl_vector* y){ //initializing 
	cinterp* s = (cinterp*)malloc(sizeof(cinterp));
	int N = x -> size;
	s -> x = gsl_vector_alloc(N);
	s -> y = gsl_vector_alloc(N);
	s -> b = gsl_vector_alloc(N);
	s -> c = gsl_vector_alloc(N-1);
	s -> d = gsl_vector_alloc(N-1);
	for(int i=0; i<N; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		gsl_vector_set(s->x, i, xi);
		gsl_vector_set(s->y, i ,yi);
	}
	gsl_vector *h = gsl_vector_alloc(N-1);
	gsl_vector *p = gsl_vector_alloc(N-1);
	for(int i=0; i<N-1; i++){
		double diffx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
		double diffy = gsl_vector_get(y, i+1)-gsl_vector_get(y,i);
		gsl_vector_set(h,i,diffx);
		gsl_vector_set(p,i,diffy/diffx);
	}
	// Tridiagonal elements
	gsl_vector *D = gsl_vector_alloc(N); // see table in notes for size
	gsl_vector *Q = gsl_vector_alloc(N-1);
	gsl_vector *B = gsl_vector_alloc(N);
	
	//Giving values to triagonal elements
	gsl_vector_set(D,0,2);
	for(int i=0; i<N-2; i++){
		double hi = gsl_vector_get(h,i);
		double hi1 = gsl_vector_get(h,i+1);
		gsl_vector_set(D,i+1,2*hi/hi1+2);
	}
	gsl_vector_set(D,N-1,2);
	
	gsl_vector_set(Q,0,1);
	for(int i=0; i<N-2; i++){
		double hi = gsl_vector_get(h,i);
		double hi1 = gsl_vector_get(h, i+1);
		gsl_vector_set(Q, i+1, hi/hi1);
	}
	
	double p0=gsl_vector_get(p,0);
	gsl_vector_set(B,0,3*p0);
	for(int i=0; i<N-2; i++){
		double pi = gsl_vector_get(p,i);
		double pi1 = gsl_vector_get(p,i+1);
		double hi = gsl_vector_get(h,i);
		double hi1 = gsl_vector_get(h,i+1);
		gsl_vector_set(B, i+1, 3*(pi+pi1*hi/hi1));
	}
	
gsl_vector_set(B,N-1, 3*gsl_vector_get(p,N-2));

	for(int i=1; i<N; i++){
	double Qi_1 = gsl_vector_get(Q, i-1);
	double Di = gsl_vector_get(D,i);
	double Di_1 = gsl_vector_get(D,i-1);
	gsl_vector_set(D,i, Di-Qi_1/Di_1);
	double Bi = gsl_vector_get(B,i);
	double Bi_1 = gsl_vector_get(B, i-1);
	gsl_vector_set(B,i, Bi-Bi_1/Di_1);
	}
	
double BN_1 = gsl_vector_get(B,N-1);
double DN_1 = gsl_vector_get(D, N-1);
gsl_vector_set(s->b, N-1, BN_1/DN_1);

	for(int i=N-2; i>=0; i--){
	double Bi = gsl_vector_get(B,i);
	double Qi = gsl_vector_get(Q,i);
	double Di = gsl_vector_get(D,i);
	double bi1 = gsl_vector_get(s->b, i+1);
	gsl_vector_set(s->b, i, (Bi-Qi*bi1)/Di);
	}	

	for(int i=0; i<N-1; i++){
	double pi = gsl_vector_get(p,i);
	double hi = gsl_vector_get(h,i);
	double bi = gsl_vector_get(s->b, i);
	double bi1 = gsl_vector_get(s->b, i+1);
	gsl_vector_set(s->c, i, (-2*bi-bi1+3*pi)/hi);
	gsl_vector_set(s->d, i, (bi+bi1-2*pi)/hi/hi);
	}
	return s;
}

	// Function that evaluates spline at given z:
	double cinterp_eval(cinterp *s, double z){
		int i = binsearch(s ->x, z);
		double h = z-gsl_vector_get(s->x,i);
		double yi = gsl_vector_get(s->y,i);
		double bi = gsl_vector_get(s->b, i);
		double ci = gsl_vector_get(s->c, i);
		double di = gsl_vector_get(s->d, i);
		double eval = yi+h*(bi+h*(ci+h*di));
	return eval;
	}

	// Function that evaluates derivative of spline at given z
	double cinterp_der(cinterp *s, double z){
		int i = binsearch(s->x, z);
		double h = z-gsl_vector_get(s->x,i);
		double bi = gsl_vector_get(s->b, i);
		double ci = gsl_vector_get(s->c, i);
		double di = gsl_vector_get(s->d, i);
		double der = bi+2*ci*h+3*di*pow(h,2);
	return der;
	}

	// Function that evaluates integral from x[0] to z
	double cinterp_integ(cinterp *s, double z){
		int j = binsearch(s->x,z);
		double area = 0;
		for(int i=0; i<j; i++){
		double h = gsl_vector_get(s->x, i+1)-gsl_vector_get(s->x,i);
		double yi = gsl_vector_get(s->y,i);
		double bi = gsl_vector_get(s->b,i);
		double ci = gsl_vector_get(s->c,i);
		double di = gsl_vector_get(s->d,i);
		area += yi*h+1./2*bi*pow(h,2)+1./3*ci*pow(h,3)+1./4*di*pow(h,4);
		}
		double hj = z-gsl_vector_get(s->x,j);
		double yj = gsl_vector_get(s->y,j);
		double bj = gsl_vector_get(s->b,j);
		double cj = gsl_vector_get(s->c,j);
		double dj = gsl_vector_get(s->d,j);
		area += yj*hj+1./2*bj*pow(hj,2)+1./3*cj*pow(hj,3)+1./4*dj*pow(hj,4);
		return area;
	}


	//Function for freeing allocated memory
	void cinterp_free(cinterp* s){
		gsl_vector_free(s->x);
		gsl_vector_free(s->y);
		gsl_vector_free(s->b);
		gsl_vector_free(s->c);
		gsl_vector_free(s->d);
		free(s);
	}

int main() {
	int N=9;
	gsl_vector *x = gsl_vector_alloc(N);
	gsl_vector *y = gsl_vector_alloc(N);
	FILE* x_file = fopen("x_points.txt", "r");
	FILE* y_file = fopen("y_points.txt", "r");
	gsl_vector_fscanf(x_file,x);
	gsl_vector_fscanf(y_file,y);


	// For GSL comparison we must have ordinary arrays:
	
	double xa[x->size];
	double ya[y->size];
	for(int i=0; i<x ->size; i++){													xa[i]=gsl_vector_get(x,i);
		ya[i]=gsl_vector_get(y,i);
	}

	cinterp* s = cinterp_alloc(x,y);
	
	// For GSL functions
	gsl_interp* cspline = gsl_interp_alloc(gsl_interp_cspline,N);
	gsl_interp_init(cspline, xa, ya, N);


	FILE* cinterp_file = fopen("cinterp.txt", "w");
	int z=0;
	double fine = 0.1;
	while(z*fine<=gsl_vector_get(x,N-1)){
		double interp_eval_gsl = gsl_interp_eval(cspline, xa, ya, z*fine, NULL);
		double interp_der_gsl = gsl_interp_eval_deriv(cspline, xa, ya, z*fine, NULL);
		double interp_integ_gsl = gsl_interp_eval_integ(cspline, xa, ya,xa[0], z*fine, NULL);
		fprintf(cinterp_file, "%10g %10g %10g %10g %10g %10g %10g\n", z*fine, cinterp_eval(s,z*fine), cinterp_der(s,z*fine), cinterp_integ(s,z*fine), interp_eval_gsl, interp_der_gsl, interp_integ_gsl);
		z++;
	}

	// freeing memory and closing files
	cinterp_free(s);
	fclose(x_file);
	fclose(y_file);
	gsl_vector_free(x);
	gsl_vector_free(y);
	fclose(cinterp_file);
	gsl_interp_free(cspline);
	return 0;
}







