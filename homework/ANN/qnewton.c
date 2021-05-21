#include<stdio.h>
#include<math.h>
#include<float.h>
static const double h = sqrt(DBL_EPSILON);
#define For(i) for(int i=0; i<d; i++)
#define ID_set(B) For(i) For(j) B[i][j]=0; For(i) B[i][i]=1;

double norm(int d, double* x){
	double sum = 0;
	For(i){
		sum += x[i]*x[i];
	}
	return sqrt(sum);
}

void gradient(
	int d,
	double f(int d, double* x),
	double* x,
	double* grad){

	For(i){
		double fx = f(d,x);
		double dx;
		double xi = x[i];
		if(fabs(xi)<h) dx = h;
		else dx = fabs(xi)*h;
		x[i]+=dx;
		grad[i]=(f(d,x)-fx)/dx;
		x[i]=xi;
	}
}

int qnewton(int d, double f(int, double*), double* x, double acc){
	int nsteps = 0;
	double B[d][d]; // inverse Hessian
	ID_set(B); // start with B being the identity
	double gradf[d];
	gradient(d, f, x, gradf); // gradient stored in gradf
	double fx = f(d,x);
	double fz;

	double maxnsteps = 10000;

	while(1){
		nsteps++;
		if(nsteps>maxnsteps) break;
		double normgradf = norm(d,gradf);
		if(normgradf<acc) break;

		double Dx[d]; // Newtons step
		For(i){ Dx[i]=0; For(j) Dx[i]-=B[i][j]*gradf[j];}
		double normx = norm(d,x);
		double normDx = norm(d,Dx);
		if(normDx<h*normx) break;

		// Start back tracking line search

		double lambda = 1;
		double z[d];
		while(1){
			For(i) z[i]=x[i]+Dx[i];
			fz = f(d,z);
			double DxTG = 0; For(i) DxTG += Dx[i]*gradf[i]; // dot product
			double alpha = 0.01;
			if(fz<fx+alpha*DxTG) break;

			if(lambda<h){ID_set(B); break;}
			lambda*=0.5;
			For(i) Dx[i]*=0.5;
		}
		
		double gradfz[d];
		gradient(d,f,z,gradfz);
		double y[d];
		For(i) y[i]=gradfz[i]-gradf[i]; // y = gradf(x+Dx)-gradf(x)
		double u[d];
		For(i){u[i]=Dx[i]; For(j) u[i]-=B[i][j]*y[j];} // u=Dx-B*y
		double uTy=0; For(i) uTy += u[i]*y[i];
		if(fabs(uTy)>1e-12){
			For(i)For(j) B[i][j] += u[i]*u[j]/uTy;
		}
		For(i) x[i]=z[i];
		For(i) gradf[i]=gradfz[i];
		fx=fz;
	}
	return nsteps;
}

