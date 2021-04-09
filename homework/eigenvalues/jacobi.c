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


// Considering only upper half part (Exercise C). 
// We write a function that do the total rotation J^T*A*J in one step, only considering the upper half of the matrix A
void JTAJ_up(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta);
	double s = sin(theta);
	int n = A -> size1;
	double app = gsl_matrix_get(A,p,p);
	double aqq = gsl_matrix_get(A,q,q);
	double apq = gsl_matrix_get(A,p,q);
	// Updating app, aqq and apq according to Eq. 10
	// By construction the angle zeroes apq after rotation
	gsl_matrix_set(A,p,q,0);
	gsl_matrix_set(A,p,p, c*c*app-2*s*c*apq+s*s*aqq);
	gsl_matrix_set(A,q,q, s*s*app+2*s*c*apq+c*c*aqq);

	// Updating the rest of upper half elements. Because of the structure of J, we only have to consider the upper half of column p and q and the right part of row p and row q
	// (see drawing on paper)
	
	for(int i=0; i<p; i++){
		double aip = gsl_matrix_get(A,i,p);
		double aiq = gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,c*aip-s*aiq); // We use that A is symmetric
		gsl_matrix_set(A,i,q,s*aip+c*aiq);
	}

	for(int i=p+1; i<q; i++){
		double api = gsl_matrix_get(A,p,i);
		double aiq = gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,p,i,c*api-s*aiq);
		gsl_matrix_set(A,i,q,s*api+c*aiq);
	}

	for(int i=q+1; i<n; i++){
		double api = gsl_matrix_get(A,p,i);
		double aqi = gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,c*api-s*aqi);
		gsl_matrix_set(A,q,i,s*api+c*aqi);
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
			// The relevant elements in the matrix are
			double apq = gsl_matrix_get(A,p,q);
			double app = gsl_matrix_get(A,p,p);
			double aqq = gsl_matrix_get(A,q,q);
			// Rotation angle
			double theta = 0.5*atan2(2*apq,aqq-app);
			double c=cos(theta);
			double s=sin(theta);
			// Calculation new elements 
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			// Check if the new element differ from the old one
			if(new_app!=app || new_aqq!=aqq)
			// If they differ we have not reached convergence and we must rotate again
			{
				changed=1;
				timesJ(A,p,q,theta);
				Jtimes(A,p,q,-theta); // Note that J^T(theta)=J(-theta)
				timesJ(V,p,q,theta);
			}
		}
		}
		}while(changed!=0);
}


// The optimized diagonalization algorithm (Exercise C):
// Same approach as above but under rotation we only change the upper half of A

void jacobi_diag_opt(gsl_matrix* A, gsl_matrix* V){
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
			JTAJ_up(A,p,q,theta);
			timesJ(V,p,q,theta);
		}
	}
	}
	}while(changed!=0);
}



