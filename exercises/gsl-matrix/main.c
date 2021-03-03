#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vector_print(char s[], gsl_vector * v) {
	printf("%s\n",s);
	for(int i=0; i < v -> size; i++)printf("%10g", gsl_vector_get(v,i));
	printf("\n");
}

int main() {
	int n = 3;
	gsl_matrix * A = gsl_matrix_alloc(n,n);
	gsl_matrix * Acopy = gsl_matrix_alloc(n,n);
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_vector * x = gsl_vector_alloc(n);
	gsl_vector * y = gsl_vector_calloc(n); // This is for the calculation of A*x (see BLAS Support)

	// inserting values in matrix A;
	double valA[3][3] = {
	{6.13 , -2.90, 5.86},
	{8.08, -6.31, -3.89},
	{-4.36, 1.00, 0.19},
	};
	for(int i=0; i < A ->size1; i++)
		for( int j=0; j < A ->size2;j++)
		{
		double Aij=valA[i][j];
		gsl_matrix_set(A,i,j,Aij);
		}
	gsl_matrix_memcpy(Acopy,A);
	
	double valb[3] = {6.23,5.37,2.29};
	for(int i=0; i< b-> size; i++)
	{
		double bi = valb[i];
		gsl_vector_set(b,i,bi);
	}
	
	gsl_linalg_HH_solve(Acopy,b,x);
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
	vector_print("x equal to:",x);
	vector_print("b equal to:",b);
	vector_print("A*x equal to:",y);

gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}


