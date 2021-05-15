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
		gsl_vector_set(grad,i,(F(x)-fx)/dx); // f(xi+dx)-f(xi)/dx
		// Get x back so its ready for the next i
		gsl_vector_set(x,i,xi);
	}
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double acc) { // Make the function an int to measure number of steps
	int n = x -> size;
	int nsteps = 0;

	gsl_matrix* B = gsl_matrix_alloc(n,n); // invers Hessian
	gsl_vector* gradF = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n); // for copying x
	gsl_vector* u = gsl_vector_alloc(n);
	gsl_vector* gradFz = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	gsl_matrix_set_identity(B); // We start with B being the identity matrix
	gradient(F,x,gradF); // calculate gradient of F, stored in gradF
	double fx = F(x); // Note that F is phi in notes

	double maxnsteps = 10000;

	while(1){ // runs until we break explicitly
		
		nsteps ++;
		if(nsteps>maxnsteps) break; 

		if(gsl_blas_dnrm2(gradF)<acc) break; // If norm of gradF is too small we stop

		gsl_blas_dgemv(CblasNoTrans, -1, B, gradF, 0, Dx); // Dx = -B*gradF Newtons step

		if(gsl_blas_dnrm2(Dx)<h*gsl_blas_dnrm2(x)) break; // step too small and we start over
		
		// Start back tracking line search
		double lambda = 1;
		while(1){ 
			gsl_vector_memcpy(z,x); // copy of x
			gsl_vector_add(z,Dx); //z=x+Dx
			double fz = F(z);
			double DxTG; 
			gsl_blas_ddot(Dx, gradF, &DxTG); // dot product of Dx and gradF
			double alpha = 0.01;
			if(fz<fx+alpha*DxTG) break; // break if Armijo condition is satisfied
			if(lambda<h){
				gsl_matrix_set_identity(B); 
				break; // if lambda is smaller than minimum we reset B to the identity and start over
			}
			lambda*=0.5;
			gsl_vector_scale(Dx,0.5);  // smaller steps
		}

		
		gradient(F,z,gradFz); // Calculate new gradient to update B
		gsl_vector_memcpy(y,gradFz); // copy
		gsl_blas_daxpy(-1, gradF, y); // y = gradF(x+Dx)-gradF(x)
		gsl_vector_memcpy(u,Dx); 
		gsl_blas_dgemv(CblasNoTrans, -1, B, y, 1, u); // u=Dx-By
		double uTy; // dot products
		gsl_blas_ddot(u,y,&uTy);
		if(fabs(uTy)>1e-12) { // SR1 update
			gsl_blas_dger(1.0/uTy,u,u,B); // B=B+uu^T/(uTy)
		}
		// Make the update
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
