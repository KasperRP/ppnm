#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>
#include<time.h>

void jacobi_diag(gsl_matrix* A, gsl_matrix* V);

void jacobi_diag_opt(gsl_matrix* A, gsl_matrix* V);

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
	


	// Now the Quantum particle in box problem
	
	FILE* quantum_eigVal_file = fopen("quantum_eigVal.txt", "w");
	FILE* quantum_eigfunc_file = fopen("quantum_eigfunc.txt","w");
	// Built Hamiltonian
	
	int n=20;
	double s=1.0/(n+1);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	for(int i=0; i<n-1; i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/s/s);

	// Diagonalize the matrix
	
	gsl_matrix* VQ = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(VQ);
	jacobi_diag(H,VQ);

	// Printing the eigenvalues of H corresponding to the energies of the system
	for (int k=0; k<n/3; k++){
		double exact = M_PI*M_PI*(k+1)*(k+1); // see paper in the map for this
		double calculated = gsl_matrix_get(H,k,k);
		fprintf(quantum_eigVal_file, "n=%i calculated=%g exact=%g\n", k+1, calculated, exact);
	}


	// Printing obtained eigenfunctions and compare to ananlytic solution sqrt(2/L)*sin(m*pi*ksi);
	double C = sqrt(2.0/(n+1));
	for(int k=0; k<3; k++){
		fprintf(quantum_eigfunc_file,"0 0 0\n");
	for(int i=0;i<n;i++){
		double ksi = (i+1.0)/(n+1);
		fprintf(quantum_eigfunc_file,"%g %g %g\n", ksi, sin((k+1)*M_PI*ksi), gsl_matrix_get(VQ, i, k)*pow(-1,k)/C);} // we multiply with (-1)^k such that the phases of our eigenfunctions agree with sin(x). And divide by C for normalization
       fprintf(quantum_eigfunc_file, "1 0 0\n");
	}       



	// For time measurements
	FILE* diag_time_file = fopen("diag_time.txt", "w");
	// Generating random  symmetric matrix with variable size
	for(int n=20; n<=200; n+=20){
		gsl_matrix* At = gsl_matrix_alloc(n,n);
		gsl_matrix* Vt = gsl_matrix_alloc(n,n);

		gsl_vector* gsl_eig = gsl_vector_alloc(n);
		gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(n);

		for(int i=0; i<n; i++){
			for(int j=i; j<n; j++){
				double Atij = (double) rand()/RAND_MAX*10;
				gsl_matrix_set(At,i,j,Atij);
				gsl_matrix_set(At,j,i,Atij);
			}
		}
		clock_t my_start = clock();
		jacobi_diag(At,Vt);
		clock_t my_stop = clock();
		double my_time = ((double) my_stop-my_start)/CLOCKS_PER_SEC;

		clock_t gsl_start = clock();
		gsl_eigen_symmv(At, gsl_eig, Vt, w);
		clock_t gsl_stop = clock();
		double gsl_time = ((double) gsl_stop-gsl_start)/CLOCKS_PER_SEC;

		clock_t my_start_opt = clock();
		jacobi_diag_opt(At,Vt);
		clock_t my_stop_opt = clock();
		double my_time_opt=((double) my_stop_opt-my_start_opt)/CLOCKS_PER_SEC;
				
		fprintf(diag_time_file, "%i %g %g %g\n", n, my_time, gsl_time, my_time_opt);

		gsl_matrix_free(At);
		gsl_matrix_free(Vt);
		gsl_vector_free(gsl_eig);
		gsl_eigen_symmv_free(w);
	}





	
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(VD);
	gsl_matrix_free(VDVT);
	gsl_matrix_free(VTV);
	gsl_matrix_free(VVT);
	
	gsl_matrix_free(H);
	gsl_matrix_free(VQ);

	fclose(jacobi_diag_file);
	fclose(quantum_eigVal_file);
	fclose(quantum_eigfunc_file);

	fclose(diag_time_file);
	return 0;
}




