#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include"binsearch.h"
#include<math.h>

// For printing results in the main function we define (see GSL Matrix exercise)

//void vector_print(char s[], gsl_vector* v) {
//printf("%s\n",s);
//for(int i=0; i < v -> size; i++) printf("%10g", gsl_vector_get(v,i));
//printf("\n");
//}

// Linear spline interpolation function
double linterp(gsl_vector* x, gsl_vector* y, double z){

int i = binsearch(x,z);
double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
double value = gsl_vector_get(y,i)+slope*(z-gsl_vector_get(x,i));
return value;
}
// Next we implement a function that calculates the integral of the linear spline from x[0] to z:

double linfunc_integ(double xi, double xf, double yi, double pi){
return yi*(xf-xi)+pi*pow(xf-xi,2)/2; // integral of straight line
}

double linterp_integ(gsl_vector* x, gsl_vector* y, double z){
int j = binsearch(x,z);
double area=0;
for(int i=0; i<j; i++){
	double slope= (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	area += linfunc_integ(gsl_vector_get(x,i), gsl_vector_get(x,i+1), gsl_vector_get(y,i), slope);
return area;
}
double slope_final= (gsl_vector_get(y,j+1)-gsl_vector_get(y,j))/(gsl_vector_get(x,j+1)-gsl_vector_get(x,j));
area += linfunc_integ(gsl_vector_get(x,j),z,gsl_vector_get(y,j),slope_final);
return area;
}



int main() {
// First the points to interpolate
int N=9; //number of tabulated points
gsl_vector* x = gsl_vector_alloc(N); 
gsl_vector* y = gsl_vector_alloc(N);
FILE* x_file = fopen("x_points.txt","r");
FILE* y_file = fopen("y_points.txt","r");
gsl_vector_fscanf(x_file,x);
gsl_vector_fscanf(y_file,y);
FILE* xy_file = fopen("xy_points.txt","w");
for(int i=0; i<=N-1; i++){
	fprintf(xy_file, "%10g %10g\n", gsl_vector_get(x,i), gsl_vector_get(y,i));
}



// Then the interpolation function


FILE* linterp_file = fopen("linterp.txt","w");
int z=0;
double fine = 0.5;
while(z*fine<=gsl_vector_get(x,N-1)){
	fprintf(linterp_file, "%10g %10g\n",z*fine, linterp(x,y,z*fine));
	z++;

}
fclose(x_file);
fclose(y_file);
fclose(xy_file);
fclose(linterp_file);
gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}



