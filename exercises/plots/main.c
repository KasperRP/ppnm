#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double myerf(double);
double mygam(double);

int main() {
	double xmin=0.1, xmax=3;
	for(double x=xmin; x<=xmax; x+=1.0/8){
		printf("%10g %10g %10g %10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),myerf(x),tgamma(x),gsl_sf_gamma(x),mygam(x));
	}
return 0;
}
