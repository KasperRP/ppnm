#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"binsearch.h"

typedef struct {gsl_vector* x, *y, *b, *c;} qinterp; // structure with relevant parameters

qinterp* qinterp_alloc(gsl_vector* x,gsl_vector* y){
	qinterp* s = (qinterp*)malloc(sizeof(qinterp)); //See the theory in the book, table 3
	int N = x -> size;
	s -> x = gsl_vector_alloc(N);
	s -> y = gsl_vector_alloc(N);
	s -> b = gsl_vector_alloc(N-1);
	s -> c = gsl_vector_alloc(N-1); // N points, N-1 intervals with polynomials
	for(int i=0; i<N; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		gsl_vector_set(s->x, i, xi);
		gsl_vector_set(s->y, i, yi); 
	}
	gsl_vector *h = gsl_vector_alloc(N-1);
	gsl_vector *p = gsl_vector_alloc(N-1);
	for(int i=0; i<N-1; i++){
		double diffx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
		double diffy = gsl_vector_get(y,i+1)-gsl_vector_get(y,i);
		gsl_vector_set(h,i,diffx);
		gsl_vector_set(p,i,diffy/diffx);
	}
	// Finding c_i by recursion
	// Recursion up
	gsl_vector_set(s->c, 0, 0); // We freely choose this (see notes)
	for(int i=0; i<N-2; i++){
		double diffp = gsl_vector_get(p,i+1)-gsl_vector_get(p,i);
		double ci=gsl_vector_get(s->c, i);
		double hi=gsl_vector_get(h,i);
		double hi1=gsl_vector_get(h,i+1);
		double ci1=(diffp-ci*hi)/hi1;
		gsl_vector_set(s->c, i+1, ci1);
	}
	// Recursion down
	gsl_vector_set(s->c, N-2, gsl_vector_get(s->c, N-2)/2);
	for(int i=N-3; i>=0; i--){
		double diffp = gsl_vector_get(p,i+1)-gsl_vector_get(p,i);
		double ci1 = gsl_vector_get(s ->c, i+1);
		double hi1 = gsl_vector_get(h,i);
		double hi = gsl_vector_get(h,i+1);
		double ci = (diffp-ci1*hi1)/hi;
		gsl_vector_set(s ->c, i, ci);
	}
	// Find b_i
	for(int i=0; i<N-1; i++){
		double pi = gsl_vector_get(p,i);
		double ci = gsl_vector_get(s->c, i);
		double hi = gsl_vector_get(h,i);
		double bi = pi-ci*hi;
		gsl_vector_set(s -> b, i, bi);
	}
	return s;
}
	// Function that evaulates interpolation in point z (see table in notes once again)
	double qinterp_eval(qinterp *s, double z){
		int i = binsearch(s->x, z);
		double h = z-gsl_vector_get(s ->x, i);
		double yi = gsl_vector_get(s ->y, i);
		double bi = gsl_vector_get(s ->b, i);
		double ci = gsl_vector_get(s ->c, i);
		double qeval = yi+h*(bi+h*ci);
		return qeval;
	}

	// Function that evaluates the derivative in point z
	double qinterp_der(qinterp *s, double z){
		int i = binsearch(s->x, z);
		double h = z-gsl_vector_get(s -> x, i);
		double bi = gsl_vector_get(s -> b, i);
		double ci = gsl_vector_get(s -> c, i);
		double qder = bi+2*ci*h;
		return qder;
	}

	// Function that evaluates integral from x[0] to z
	double qinterp_integ(qinterp *s, double z){
		int j = binsearch(s->x, z);
		double area = 0;
		for(int i=0; i<j; i++){ // summing over all intervals except the last one
		double h = gsl_vector_get(s ->x, i+1)-gsl_vector_get(s->x, i);
		double yi = gsl_vector_get(s ->y, i);
		double bi = gsl_vector_get(s ->b, i);
		double ci = gsl_vector_get(s ->c, i);
		area += yi*h+1./2*bi*pow(h,2)+1./3*ci*pow(h,3);
		}	
		double hj = z-gsl_vector_get(s->x, j); // and then the last one
		double yj = gsl_vector_get(s ->y, j);
		double bj = gsl_vector_get(s ->b, j);
		double cj = gsl_vector_get(s ->c, j);
		area += yj*hj+1./2*bj*pow(hj,2)+1./3*cj*pow(hj,3);
		return area;
	}

	//Function for freeing allocated memory
	void qinterp_free(qinterp* s){
		gsl_vector_free(s ->x);
		gsl_vector_free(s ->y);
		gsl_vector_free(s ->b);
		gsl_vector_free(s ->c);
		free(s);
	}


	int main() {
	int N=9;
	gsl_vector *x = gsl_vector_alloc(N);
	gsl_vector *y = gsl_vector_alloc(N);
	FILE* x_file = fopen("x_points.txt","r");
	FILE* y_file = fopen("y_points.txt", "r");
	gsl_vector_fscanf(x_file,x);
	gsl_vector_fscanf(y_file,y);
	
	qinterp* s = qinterp_alloc(x,y); 


	FILE* qinterp_file = fopen("qinterp.txt", "w");
	int z=0;
	double fine = 0.1;
	while(z*fine<=gsl_vector_get(x,N-1)){
		fprintf(qinterp_file, "%10g %10g %10g %10g\n", z*fine, qinterp_eval(s, z*fine), qinterp_der(s, z*fine), qinterp_integ(s,z*fine));
		z++;
	}
	
	// freeing memory and closing files using the above defined function
	qinterp_free(s);
	fclose(x_file);
	fclose(y_file);
	gsl_vector_free(x);
	gsl_vector_free(y);
	fclose(qinterp_file);
	return 0;
}





