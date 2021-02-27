#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f(double t, void* params){
	double z = *(double*)params;
	double f = 2/sqrt(M_PI)*exp(-z*t*t);
	return f;
}

double myerf(double x){
	double z=1;
	gsl_function F;
	F.function = &f;
	F.params = (void*) &z;
	int limit = 999;
	gsl_integration_workspace *w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0, acc=1e-6, eps=1e-6,result,error;
	gsl_integration_qags(&F,a,x,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main() {
	double x;
	double xmin=-2, xmax=2;
	for(x=xmin;x<=xmax;x+=1.0/8)
		printf("%10g %10g\n",x,myerf(x));
return 0;
}	
