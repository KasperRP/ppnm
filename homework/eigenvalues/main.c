#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void jacobi_diag(gsl_matrix* A, gsl_matrix* V);

void matrix_print(gsl_matrix* A, FILE* list){
	int m = A -> size1;
	int n = A -> size2;
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			double Aij = gsl_matrix_get(A,i,j);
			fprintf(list, "%.2f  ", Aij);
		}
		fprintf(list, "\n");
	}
}


int main(){
	
	//Generating file for checking jacobi_diag function
	FILE* jacobi_diag_file = fopen("jacobi_diag_test.txt", "w");

	// Generating random real symmetric matrix A
	int N = 5;
	gsl_matrix* A = gsl_matrix_alloc(N,N);
	for(int i=0; i<N; i++){
		for(int j=i; j<N; j++){
			double Aij = (double) rand()/RAND_MAX*10;
			gsl_matrix_set(A,i,j,Aij);
			gsl_matrix_set(A,j,i,Aij);
		}
	}

	// printing the relevant matrices
	fprintf(jacobi_diag_file, "Random real symmetric matrix A:\n");
	matrix_print(A, jacobi_diag_file);

	gsl_matrix* V = gsl_matrix_alloc(N,N);
	gsl_matrix_set_identity(V); // See page 36 in notes
	jacobi_diag(A, V); // Stores the diagonalmatrix in A
	fprintf(jacobi_diag_file, "\nDiagonal matrix containing eigenvalues of A:\n");
	matrix_print(A, jacobi_diag_file);
	
	
	
	// Now I check that VDV^T=A
	
	gsl_matrix* VD = gsl_matrix_alloc(N,N);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, A, 0, VD);

	gsl_matrix* VDVT = gsl_matrix_alloc(N,N);
	
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, VD, V, 0, VDVT);

	fprintf(jacobi_diag_file, "\nThe product VDV^T. Should be equal to original A:\n");
	matrix_print(VDVT, jacobi_diag_file);

	gsl_matrix* VTV = gsl_matrix_alloc(N,N);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VTV);
	gsl_matrix* VVT = gsl_matrix_alloc(N,N);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, V, V, 0, VVT);
	
	fprintf(jacobi_diag_file, "\nHere comes the products V^T*V and VV^T. Both should be the identity:\n");
	matrix_print(VTV, jacobi_diag_file);
	fprintf(jacobi_diag_file, "\n \n");
	matrix_print(VVT, jacobi_diag_file);


	
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(VD);
	gsl_matrix_free(VDVT);
	gsl_matrix_free(VTV);
	gsl_matrix_free(VVT);
	fclose(jacobi_diag_file);
	return 0;
}




