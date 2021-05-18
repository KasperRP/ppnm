#include<stdio.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

int qnewton(double F(gsl_vector*), gsl_vector* x, double acc); // minimization function


