#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// Calling fitting function

void lsfunc(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double funs(int, double), int m, gsl_vector* c, gsl_vector* dc, gsl_matrix* Cov);

// Defining function for linear fitting

double linfunc(int i, double x){
	if(i==0) return 1;
	else if(i==1) return x;
	else return NAN;
}



int main(){
	
	int N = 9; // number of data points
	gsl_vector* t = gsl_vector_alloc(N);
	gsl_vector* y = gsl_vector_alloc(N);
	FILE* t_file = fopen("time.txt", "r");
	FILE* act_file = fopen("activity.txt", "r");
	gsl_vector_fscanf(t_file, t);
	gsl_vector_fscanf(act_file, y);

	gsl_vector* logy = gsl_vector_alloc(N);
	gsl_vector* d_logy = gsl_vector_alloc(N);
	for(int i=0; i<N; i++){
		double yi = gsl_vector_get(y,i);
		double dyi = yi/20; // Uncertainty in y as given in exercise
		double log_yi = log(yi);
		double d_log_yi = dyi/yi;
		gsl_vector_set(logy, i , log_yi); // We are going to plot ln(y) and its uncertainty
		gsl_vector_set(d_logy, i, d_log_yi); 
	}
	
	// Making file with datapoints and errors
	FILE* data_file = fopen("data.txt", "w");
	
	for(int i=0; i<N; i++){
		fprintf(data_file, "%g %g %g\n", gsl_vector_get(t,i), gsl_vector_get(logy, i), gsl_vector_get(d_logy,i));
	}



	
	int M = 2; // Corresponding to linear fitting
	gsl_vector* c = gsl_vector_alloc(M);
	gsl_vector* dc = gsl_vector_alloc(M);
	gsl_matrix* Cov = gsl_matrix_alloc(M,M);

	//Now we do the fitting and estimate the half-lifetime
	FILE* half_life_time = fopen("half_life_time.txt", "w");

	lsfunc(t, logy, d_logy, linfunc, 2, c, dc, Cov);
	double c0 = gsl_vector_get(c,0);
	double c1 = gsl_vector_get(c,1);
	double dc0 = gsl_vector_get(dc,0);
	double dc1 = gsl_vector_get(dc,1);
	double k = -c1;
	double half_lf = log(2)/k;
	double err_k=dc1;
	double err_half_lf = log(2)/(k*k)*err_k; // using error-propagation law

	fprintf(half_life_time, "The coefficients in my linear model are found to be c(0)=%g and c(1)=%g with uncertainties dc(0)=%g and dc(1)=%g.\nFrom this I find the decay rate to be k=-c(1)=%g per day with uncertaintity err_k=%g per day.\nThe estimated half-lifetime is then T=log(2)/k=%g days with uncertainty (using error propagation law) log(2)/k^2*err_k=%g days.\nThe modern value of the half-lifetime of Ra-224 is 3.632 days according to Wikipedia, so even my lowest possible value is above the modern one.\n", c0, c1, dc0, dc1 , k, err_k, half_lf, err_half_lf);
	
	FILE* fit = fopen("fit.txt", "w");
	int z=0;
	double fine=0.1;
	while(z*fine <=20){
		fprintf(fit, "%g %g %g %g\n", z*fine, c0+c1*z*fine, (c0-dc0)+(c1-dc1)*z*fine, (c0+dc0)+(c1+dc1)*z*fine);
		z++;
	}
	


	gsl_vector_free(t);
	gsl_vector_free(y);
	gsl_vector_free(logy);
	gsl_vector_free(d_logy);
	fclose(t_file);
	fclose(act_file);
	fclose(data_file);


	gsl_vector_free(c);
	gsl_vector_free(dc);
	gsl_matrix_free(Cov);
	fclose(half_life_time);
	fclose(fit);
	return 0;
}
