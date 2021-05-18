#include<stdio.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void gradient(double F(gsl_vector*), gsl_vector* x, gsl_vector* grad);
int qnewton(double F(gsl_vector*), gsl_vector* x, double acc);

int amoeba(int dim, double f(double*), double** simplex, double simplex_size_goal); // For part C

double Rosenbrock(gsl_vector* r){  // (A)
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}
double RosenbrockC(double* r){ // (C) (No use of GSL vectors)
	return (1-r[0])*(1-r[0])+100*(r[1]-r[0]*r[0])*(r[1]-r[0]*r[0]);
}

double Himmelblau(gsl_vector* r){ // (A)
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

double HimmelblauC(double* r){
	return (r[0]*r[0]+r[1]-11)*(r[0]*r[0]+r[1]-11)+(r[0]+r[1]*r[1]-7)*(r[0]+r[1]*r[1]-7);
}

// Breit-Wigner function(B)
double BW(double Energy, gsl_vector* p){ 
	double m = gsl_vector_get(p, 0);
	double Gamma = gsl_vector_get(p, 1);
	double A = gsl_vector_get(p, 2);
	return A/((Energy-m)*(Energy-m)+Gamma*Gamma/4);
}

// Deviation function (B)
int ndata;
gsl_vector* E;
gsl_vector* sigma;
gsl_vector* dsigma;

double D(gsl_vector* p){
	double sum =0;
	for(int i=0; i<ndata; i++){
	double Ei = gsl_vector_get(E,i);
	double sigmai = gsl_vector_get(sigma, i);
	double dsigmai = gsl_vector_get(dsigma,i);
	double Di = (BW(Ei, p)-sigmai)*(BW(Ei,p)-sigmai)/dsigmai/dsigmai;
	sum += Di;
	}
	return sum;
}
	



int main(){

	{ // A 
	FILE* Exc_A = fopen("Exc_A.txt", "w");
	// Minimum of Rosenbrock
	int n = 2;
	gsl_vector* x0 = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector_set(x0, 0, -2);
	gsl_vector_set(x0, 1, 8); // initial guess
	double acc = 1e-4;
	gsl_vector_memcpy(x,x0);
	int steps = qnewton(Rosenbrock, x, acc);
	fprintf(Exc_A, "Rosenbrock's valley function:\n");
	fprintf(Exc_A, "Tolerance = %g\n", acc);
	fprintf(Exc_A, "Initial guess: (%g ,%g)\n", gsl_vector_get(x0,0),gsl_vector_get(x0,1) );
	fprintf(Exc_A, "Found minimum: (%g,%g)\n", gsl_vector_get(x,0), gsl_vector_get(x,1) );
	fprintf(Exc_A, "Value at min: %g\n", Rosenbrock(x));
	fprintf(Exc_A, "Number of steps : %i\n", steps);
	
	// Himmelblau
	gsl_vector* y0 = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector_set(y0, 0, 4);
	gsl_vector_set(y0, 1, 2); // initial guess
	double eps = 0.00001;
	gsl_vector_memcpy(y,y0);
	int stepss = qnewton(Himmelblau, y, eps);
	fprintf(Exc_A, "\nHimmelblau's function: \n");
	fprintf(Exc_A, "Tolerance = %g\n", eps);
	fprintf(Exc_A, "Initial guess: (%g,%g)\n", gsl_vector_get(y0,0), gsl_vector_get(y0,1) );
	fprintf(Exc_A, "Found minimum: (%g,%g)\n", gsl_vector_get(y,0), gsl_vector_get(y,1) );
	fprintf(Exc_A, "Value at min: %g\n", Himmelblau(y));
	fprintf(Exc_A, "Number of steps: %i\n", stepss);
	
	fclose(Exc_A);
	gsl_vector_free(x0);
	gsl_vector_free(x);
	gsl_vector_free(y0);
	gsl_vector_free(y);
	}

	{ // Part B
	FILE* Exc_B = fopen("Exc_B.txt", "w");
	// Getting Higgs data
	ndata = 30;
	E = gsl_vector_alloc(ndata);
	sigma = gsl_vector_alloc(ndata);
	dsigma = gsl_vector_alloc(ndata);
	FILE* Energy_file = fopen("Energy.txt", "r");
	FILE* sigma_file = fopen("sigma.txt", "r");
	FILE* dsigma_file = fopen("dsigma.txt", "r");
	gsl_vector_fscanf(Energy_file, E);
	gsl_vector_fscanf(sigma_file, sigma);
	gsl_vector_fscanf(dsigma_file, dsigma);
	fclose(Energy_file);
	fclose(sigma_file);
	fclose(dsigma_file);
	



	int N = 3; // Number of parameters (m, gamma, A)
	gsl_vector* param0 = gsl_vector_alloc(N);
	gsl_vector* param = gsl_vector_alloc(N);	
	gsl_vector_set(param0,0,125);
	gsl_vector_set(param0,1,1);
	gsl_vector_set(param0,2,1); // initial guess

	gsl_vector_memcpy(param, param0);
	double tol = 1e-3;
	int nsteps = qnewton(D, param, tol);
	fprintf(Exc_B, "Higgs Boson: \n");
	fprintf(Exc_B, "Tolerance = %g\n", tol);
	fprintf(Exc_B, "Initial guess: (m, Gamma, A) = (%g,%g,%g)\n", gsl_vector_get(param0,0), gsl_vector_get(param0,1), gsl_vector_get(param0,2) );
	fprintf(Exc_B, "Estimate values: (m, Gamma, A) =(%g,%g,%g)\n", gsl_vector_get(param,0), gsl_vector_get(param,1), gsl_vector_get(param,2) );
	fprintf(Exc_B, "Number of steps: %i\n", nsteps);

	fclose(Exc_B);
	gsl_vector_free(param0);
	gsl_vector_free(param);
	}

	{// Part C
	FILE* Exc_C = fopen("Exc_C.txt", "w");
	//Minimum of Rosenbrock's valley function
	int dim = 2;
	double* simplex[] = { // initial value of simplex. How hard it is for convergence depends on this value

	(double[]) {1.0, 2.0},
	(double[]) {0.0, 2.0},
	(double[]) {8.0, 0.0}
	};
	double size_goal = 1e-4;
	//int hi = 0;
	int lo = 0;

	int steps = amoeba(dim, RosenbrockC, simplex, size_goal);
	fprintf(Exc_C, "Rosenbrock's valley function: \n");
	fprintf(Exc_C, "Simplex size goal = %g\n", size_goal);
	fprintf(Exc_C, "Found minimum at (%g, %g)\n", simplex[lo][0], simplex[lo][1]);
	fprintf(Exc_C, "Number of steps : %i\n", steps);
	
	double* Simplex[] = {
	
	(double[]) {1.0, 2.0},
	(double []) {0.0, 2.0},
	(double []) {8.0, 0.0}
	};
	
	

	int stepss = amoeba(dim, HimmelblauC, Simplex, size_goal);
	fprintf(Exc_C, "Himmelblau function: \n");
	fprintf(Exc_C, "Simplex size goal = %g\n", size_goal);
	fprintf(Exc_C, "Found minimum at (%g,%g)\n", Simplex[lo][0], Simplex[lo][1]);
	fprintf(Exc_C, "Number of steps: %i\n", stepss);

	fclose(Exc_C);
	}
		

	


	return 0;
}

