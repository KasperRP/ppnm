#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

// Function that multiply a given matrix A with J from right

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta);
	double s=sin(theta);
	int n = A -> size1;
	for(int i=0; i<n; i++){
		double aip = gsl_matrix_get(A,i,p);
		double aiq = gsl_matrix_get(A,i,q);
		double new_aip = c*aip-s*aiq; 
		double new_aiq = s*aip+c*aiq;
		gsl_matrix_set(A, i, p, new_aip);
		gsl_matrix_set(A, i, q, new_aiq);
	}
}

// Function that multiply A with J from left

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta);
	double s=sin(theta);
	int m = A -> size2;
	for(int j=0; j<m; j++){
		double apj = gsl_matrix_get(A,p,j);
		double aqj = gsl_matrix_get(A,q,j);
		double new_apj = c*apj+s*aqj;
		double new_aqj = -s*apj+c*aqj;
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
	}
}

// Jacobi eigenvalue algorithm for matrix diagonalization

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	int n = A -> size1;
	int changed;
	do{
		changed=0;
		for(int p=0; p<n-1; p++){
			for(int q=p+1; q<n; q++){
			double apq = gsl_matrix_get(A,p,q);
			double app = gsl_matrix_get(A,p,p);
			double aqq = gsl_matrix_get(A,q,q);
			double theta = 0.5*atan2(2*apq,aqq-app);
			double c=cos(theta);
			double s=sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			if(new_app!=app || new_aqq!=aqq)
			{
				changed=1;
				timesJ(A,p,q,theta);
				Jtimes(A,p,q,-theta);
				timesJ(V,p,q,theta);
			}
		}
		}
		}while(changed!=0);
}










