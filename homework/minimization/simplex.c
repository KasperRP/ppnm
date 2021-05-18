#include<stdio.h>
#include<math.h>

void contraction(
		double* highest, double* centroid, double* contracted, int dim){
	for(int k=0; k<dim; k++) contracted[k]=0.5*centroid[k]+0.5*highest[k];
}

void reflection(
		double* highest, double* centroid, double* reflected, int dim){
	for(int k=0; k<dim; k++) reflected[k]=2*centroid[k]-highest[k];
}

void expansion(
		double* highest, double* centroid, double* expanded, int dim){
	for(int k=0; k<dim; k++) expanded[k]=3*centroid[k]-2*highest[k];
}

void reduction( double** simplex, int dim, int lo){
	for(int k=0; k<dim+1; k++) if(k!=lo) for(int i=0; i<dim; i++)
		simplex[k][i]=0.5*(simplex[k][i]+simplex[lo][i]);
}

double distance(double* a, double* b, int dim){
	double s=0;
	for(int k=0; k<dim; k++) s += (b[k]-a[k])*(b[k]-a[k]);
	return sqrt(s);
}
// size of simplex
double size(double** simplex, int dim){
	double s=0;
	for(int j=1; j<dim+1; j++){
		double d = distance(simplex[0], simplex[j], dim);
		if(d>s) s=d; }
	return s;
}

void simplex_update(
		int dim,
		double** simplex,
		double* f_val,
		int* hi,
		int* lo,
		double* centroid){

	*hi=0;
	*lo=0;
	double fhi = f_val[0];
	double flo = f_val[0];

	// Find highest and lowest value
	for(int k=1; k<dim+1; k++){
		double fx = f_val[k];
		if(fx>fhi) {fhi=fx; *hi=k;}
		if(fx<flo) {flo=fx;*lo=k;}
	}

	//Finding centroid
	for(int i=0; i<dim; i++){
		double sum=0;
		for(int j=0; j<dim+1;j++) if(j != *hi) sum += simplex[j][i];
			centroid[i]=sum/dim;
		
	}
}
	
// Initiate
void simplex_init(
		int dim,
		double fun(double*),
		double** simplex,
		double* f_val,
		int* hi,
		int* lo,
		double* centroid){
	// function value at vertices
	for(int j=0; j<dim+1; j++) f_val[j]=fun(simplex[j]);
	simplex_update(dim, simplex, f_val, hi, lo, centroid);
}


int amoeba(
	int dim,
	double f(double*),
	double** simplex,
	double simplex_size_goal){
	
	int nsteps = 0;

	int hi, lo;
	double centroid[dim];
	double f_val[dim+1];
	double p1[dim], p2[dim];
	
	// Start by initiating simplex
	simplex_init(dim, f, simplex, f_val, &hi, &lo, centroid);

	while(1){ // loop until we break explicitly
		
		if(size(simplex,dim)<simplex_size_goal) break;
		
		// Get ready for reflection
		simplex_update(dim, simplex, f_val, &hi, &lo, centroid);
		reflection(simplex[hi], centroid, p1, dim);
		double f_re = f(p1);

		if(f_re<f_val[lo]){ //if reflection OK, try expansion
		expansion(simplex[hi], centroid, p2, dim);
		double f_ex=f(p2);

		if(f_ex<f_re) { // then we accept expansion
		for(int i=0; i<dim; i++) simplex[hi][i]=p2[i]; f_val[hi]=f_ex;
	       	}

		else{ // we reject expansion and accept reflection

			for(int i=0; i<dim; i++) simplex[hi][i]=p1[i]; f_val[hi]=f_re;}
		}

		// if reflection is not good, we try the following
		else{ 
			if(f_re<f_val[hi]){ // then we accept reflection anyway
				for(int i=0; i<dim; i++) simplex[hi][i]=p1[i]; f_val[hi]=f_re;
			}

			else{ // else we try contraction
				contraction(simplex[hi], centroid, p1, dim); 
				double f_co = f(p1);

				if(f_co<f_val[hi]){ // then we accept contraction
					for(int i=0; i<dim; i++) simplex[hi][i]=p1[i]; f_val[hi]=f_re;
				}
				else{ // try reduction
					reduction(simplex, dim, lo);

					simplex_init(dim, f, simplex, f_val, &hi, &lo, centroid);}}}
	
		nsteps++;
		}

		return nsteps;
}


		





	
