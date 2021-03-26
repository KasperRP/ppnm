#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// Defining dot product and norm to ease the notation

double dot(gsl_vector* x, gsl_vector* y){
	double x_dot_y;
	gsl_blas_ddot(x, y, &x_dot_y); //GSL function that calculates dot product
	return x_dot_y;
}

double norm(gsl_vector* x){
	return sqrt(dot(x,x));
}
 // Decomposition using modified Gram-Schmidt algorithm. It takes a matrix A and decompose it A=Q*R, where Q 
 // (here stored as A) is orthogonal and R is upper triangular)
 
void GS_decomp(gsl_matrix* A, gsl_matrix* R){
	int n = A->size2;
	for(int i=0; i<n; i++){
		gsl_vector_view ai_view = gsl_matrix_column(A,i); // view of column i of A as vector
		gsl_vector* ai = &ai_view.vector;
		double Rii = norm(ai);
		gsl_matrix_set(R,i,i,Rii); // Setting the diagonalelements of R
		gsl_vector_scale(ai, 1/Rii); // normalize ai

		for(int j=i+1; j<n; j++){
			gsl_vector_view aj_view = gsl_matrix_column(A,j);
			gsl_vector* aj = &aj_view.vector;
			double Rij = dot(ai, aj);
			gsl_matrix_set(R,i,j,Rij);
			gsl_blas_daxpy(-Rij,ai, aj); // aj -> aj-ai*Rij
		}
	}
}



// Implementing in-place substitution (see notes)
// It takes the upper triangular matrix R and a vector (Q^T*b, here denoted x), solves the equation R*x=Q^T*b,
// and puts the solution into the vector x

void backsub(gsl_matrix* R, gsl_vector* x){
	int n = x -> size;
	for(int i=n-1; i>=0; i--){
		double xi = gsl_vector_get(x,i);
		for(int k=i+1; k<n; k++){
			xi-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);}
		gsl_vector_set(x,i,xi/gsl_matrix_get(R,i,i));
	}
}



void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
// We will solve QRx=b. First we need the equation to have the form Rx=Q^T*b
gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x); // Calculates Q^T*b, contains the result in x
backsub(R,x); // solves R*x= Q^T*b (initially x=Q^T*b, but replaced by the solution x after running)
}


void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* X){
// Finding inverse of A=QR by solving QR*xi=ei. Inverse is stored in X
int n=Q -> size2;
gsl_vector* ei = gsl_vector_alloc(n);
for(int i=0; i<n; i++){
	gsl_vector_set_basis(ei,i); //Making ei the i'th unit vector
	gsl_vector_view xi_view = gsl_matrix_column(X,i);
	gsl_vector* xi = &xi_view.vector;
	GS_solve(Q, R, ei, xi); //Solves QR*xi=ei. Result ends in X because xi is defined as column of X
}
	gsl_vector_free(ei);
	}

