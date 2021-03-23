#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// Calling functions from gs.c

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

// Defining functions for printing matrices and vectors

void matrix_print(gsl_matrix* A, FILE* list){
	int m = A -> size1;
	int n = A -> size2;
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			double Aij = gsl_matrix_get(A,i,j);
			fprintf(list,"%g  ", Aij);
		}
		fprintf(list, "\n");
	}
}

void vector_print(gsl_vector* b, FILE* list){
	int n = b -> size;
	for(int i=0; i<n; i++){
		double bi = gsl_vector_get(b,i);
		fprintf(list,"%g\n", bi);
	}
}

int main(){
	
	//Generate file for checking GS_decomp function
	FILE* GS_decomp_file = fopen("GS_decomp_test.txt", "w");

	// Generates random matrix
	int m = 5;
	int n = 3;
	gsl_matrix* A = gsl_matrix_alloc(m,n);
	for(int i=0; i<m; i++){
		for(int j=0; j<n;j++){
			double Aij = (double) rand()/RAND_MAX*10; // divide by RAND_MAX to avoid insanely large numbers
			gsl_matrix_set(A,i,j,Aij);
		}
	}
	fprintf(GS_decomp_file, "Random generates matrix A:\n");
	matrix_print(A, GS_decomp_file);


	// Factorize A=QR;
	gsl_matrix* R = gsl_matrix_alloc(n,n);

	GS_decomp(A,R);
	fprintf(GS_decomp_file, "\nHere we have Q:\n");
	matrix_print(A, GS_decomp_file);
	fprintf(GS_decomp_file, "\nAnd R. It should be upper triangular:\n");
	matrix_print(R, GS_decomp_file);

	gsl_matrix* QTQ = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QTQ); // Calculate Q^T*Q
	fprintf(GS_decomp_file, "\nHere comes Q^T*Q. It should be the identity\n");
	matrix_print(QTQ, GS_decomp_file);

	gsl_matrix* QR = gsl_matrix_alloc(m,n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
	fprintf(GS_decomp_file, "\nQ*R. This should be equal to A:\n");
	matrix_print(QR, GS_decomp_file);



	//Generate file for checking GS_solve
	FILE* GS_solve_file = fopen("GS_solve_test.txt", "w");
	// Generate random square matrix and vector b
	int N=4;
	gsl_matrix* A2 = gsl_matrix_alloc(N,N);
	gsl_vector* b = gsl_vector_alloc(N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			double A2ij=(double) rand()/RAND_MAX*10;
			gsl_matrix_set(A2,i,j,A2ij);
		}
	}
	for(int i=0; i<N; i++){
		double bi = (double) rand()/RAND_MAX*10;
		gsl_vector_set(b,i,bi);
	}

	fprintf(GS_solve_file, "Random generated matrix A:\n");
	matrix_print(A2, GS_solve_file);

	// Factorize A=QR:
	gsl_matrix* R2 = gsl_matrix_alloc(N,N);
	GS_decomp(A2,R2);
	fprintf(GS_solve_file, "\nHere we have Q:\n");
	matrix_print(A2, GS_solve_file);
	fprintf(GS_solve_file, "\nAnd R, the upper triangular one:\n");
	matrix_print(R, GS_solve_file);

	//Solve QRx=b
	fprintf(GS_solve_file, "\nHere is the random generated vector b: \n");
	vector_print(b, GS_solve_file);
	gsl_vector* x = gsl_vector_alloc(N);
	GS_solve(A2, R2, b, x);
	fprintf(GS_solve_file, "\nHere is the solution x to QRx=b:\n");
	vector_print(x, GS_solve_file);
	gsl_matrix* QR2 = gsl_matrix_alloc(N,N);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A2, R2, 0, QR2); // Calculates matrix product Q*R
	gsl_vector* QRx = gsl_vector_alloc(N);
	gsl_blas_dgemv(CblasNoTrans, 1, QR2, x, 0, QRx); // Calculates Q*R*x
	fprintf(GS_solve_file, "\nHere is QRx. It should be equal to b:\n");
	vector_print(QRx, GS_solve_file);

	



	
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(QTQ);
	gsl_matrix_free(QR);
	fclose(GS_decomp_file);	
	

	gsl_matrix_free(A2);
	gsl_matrix_free(R2);
	gsl_matrix_free(QR2);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(QRx);
	fclose(GS_solve_file);
	return 0;
}






	
