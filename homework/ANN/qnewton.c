#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
static const double h = sqrt(DBL_EPSILON);

void gradient(
		double F(gsl_vector*),
		gsl_vector* x,
		gsl_vector* grad){
	double m = x -> size;
	for(int i=0; i<m; i++){
		double fx = F(x);
		double dx;
		double xi = gsl_vector_get(x,i);
		if(fabs(xi)<h) dx=h;
		else dx=fabs(xi)*h;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double acc) {
	int n = x->size;
	int nsteps =0;

	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_vector* gradF = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);
	gsl_vector* gradFz = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	gsl_matrix_set_identity(B);
	gradient(F,x,gradF);
	double fx = F(x);

	double maxnsteps = 10000;

	while(1){

		nsteps ++;
		if(nsteps>maxnsteps) break;

		if(gsl_blas_dnrm2(gradF)<acc) break;

		gsl_blas_dgemv(CblasNoTransm -1, B, gradF, 0, Dx);

		if(gsl_blas_dnrm2(Dx)<h*gsl_blas_dnrm2(x)) break;

		
		double lambda = 1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			double fz = F(z);
			double DxTG;
			gsl_blas_ddot(Dx, gradF, &DxTG);
			double alpha = 0.01;
			if(fz<fx+alpha*DxTG) break;
			if(lambda<h){
				gsl_matrix_set_identity(B);
				break;
			}
			lambda*=0.5;
			gsl_vector_scale(Dx,0.5);
		}

		gradient(F,z,gradFz);
		gsl_vector_memcpy(y,gradFz);
		gsl_blas_daxpy(-1, gradF, y);
		gsl_vector_memcpy(u,Dx);
		gsl_blas_dgemv(CblasNoTrans, -1, B, y, 1, u);
		double uTy;
		gsl_blas_ddot(u,y,&uTy);
		if(fabs(uTy)>1e-12) {
			gsl_blas_dger(1.0/uTy,u,u,B);
		}

		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gradF,gradFz);
		fx=F(x);}

		gsl_matrix_free(B);
		gsl_vector_free(gradF);
		gsl_vector_free(Dx);
		gsl_vector_free(z);
		gsl_vector_free(u);
		gsl_vector_free(gradFz);
		gsl_vector_free(y);

		return nsteps;
}

