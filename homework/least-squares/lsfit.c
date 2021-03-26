#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// Calling functions from gs.c:

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* X);

// Function returning the coefficients in linear combination of fitting functions

void lsfunc(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double funs(int, double), int m, gsl_vector* c, gsl_vector* dc, gsl_matrix* Cov){  // m is number of fitting functions

	//defining matrix A and vector b according to (7) in notes:
	int n = x -> size;

	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i=0; i<n; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		double dyi = gsl_vector_get(dy,i);
		double bi = yi/dyi;
		gsl_vector_set(b,i,bi);
		for(int k=0; k<m; k++){
			double Aik = funs(k,xi)/dyi;
			gsl_matrix_set(A,i,k,Aik);
		}
	}
	// decompose A into A=QR:
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	GS_decomp(A,R);
	// Then we can solve QRc=b using the function GS_solve
	GS_solve(A, R, b, c); // stores solution in vector c
	
// Now we find the covariance matrix given by (15) in notes. GS_inverse(Q, R, X) stores the inverse of Q*R in X. We are interested in the inverse of R so we let Q=I:
	gsl_matrix* I = gsl_matrix_alloc(m,m);
	gsl_matrix_set_identity(I);
	gsl_matrix* R_inv = gsl_matrix_alloc(m,m);
	GS_inverse(I, R, R_inv); //calculates inverse of R and stores it in Riv

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, R_inv, R_inv, 0, Cov); // Cov = R^-1*(R^-1)^T

	// Finally we find the uncertainties of c given by Eq. (16)
	for(int k=0; k<m; k++){
	double Covkk = gsl_matrix_get(Cov, k, k);
	double dck = sqrt(Covkk);
	gsl_vector_set(dc, k, dck);
	}

	
	


	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(b);

	gsl_matrix_free(I);
	gsl_matrix_free(R_inv);
	


}
	
