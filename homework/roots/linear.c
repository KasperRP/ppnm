#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// Linear equation solver from exercise "Linear Equations"

double dot(gsl_vector* x, gsl_vector* y){
	double x_dot_y;
	gsl_blas_ddot(x,y,&x_dot_y);
	return x_dot_y;
}

double norm(gsl_vector* x){
	return sqrt(dot(x,x));
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R){
	int n = A->size2;
	for(int i=0; i<n; i++){
	gsl_vector_view ai_view = gsl_matrix_column(A,i);
	gsl_vector* ai = &ai_view.vector;
	double Rii = norm(ai);
	gsl_matrix_set(R,i,i,Rii);
	gsl_vector_scale(ai, 1/Rii);

	for(int j=i+1; j<n; j++){
		gsl_vector_view aj_view = gsl_matrix_column(A,j);
		gsl_vector* aj = &aj_view.vector;
		double Rij = dot(ai, aj);
		gsl_matrix_set(R,i,j,Rij);
		gsl_blas_daxpy(-Rij,ai,aj);
	}
	}
}

void backsub(gsl_matrix* R, gsl_vector* x){
	int n = x -> size;
	for(int i=n-1; i>=0; i--){
		double xi = gsl_vector_get(x,i);
		for(int k=i+1; k<n; k++){
			xi-= gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);}
		gsl_vector_set(x,i,xi/gsl_matrix_get(R,i,i));
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
	backsub(R,x);
}
