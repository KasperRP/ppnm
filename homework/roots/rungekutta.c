#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void rkstep23(
	void f(int n, double x, double* y, double* dydx),
	int n,
	double x,
	double* yx,
	double h,
	double* yh,
	double* dy) {

	double k1[n], k2[n], k3[n], k4[n];
	double yt[n];

	f(n,x,yx,k1);
	for(int i=0; i<n; i++) {yt[i]=yx[i]+1./2*k1[i]*h;}

	f(n, x+1./2*h, yt, k2);
	for(int i=0; i<n; i++) {yt[i]=yx[i]+3./4*k2[i]*h;}

	f(n, x+3./4*h, yt, k3);
	for(int i=0; i<n; i++) {yh[i]=yx[i]+(2./9*k1[i]+1./3*k2[i]+4./9*k3[i])*h;}

	f(n, x+h, yh, k4);
	for(int i=0; i<n; i++){
		yt[i]=yx[i]+(7./24*k1[i]+1./4*k2[i]+1./3*k3[i]+1./8*k4[i])*h;
		dy[i]=yt[i]-yh[i];
	}
}

void driver(
		void f(int n, double x, double* y, double* dydx),
		int n,
		double a,
		double* ya,
		double b,
		double* yb,
		double h,
		double acc,
		double eps,
		char* path
	   )
{
	FILE* list = fopen(path, "w");
	double x;
	double y[n];
	double sum;
	double err;
	double normy;
	double tol;
	double yh[n];
	double dy[n];

	x=a;
	fprintf(list, "%20g", x);
	for(int i=0; i<n; i++){
		y[i]=ya[i];
		fprintf(list, "%20g", y[i]);
	}
	fprintf(list, "\n");

	while(x<b){

		if(x+h>b) h=b-x;

	rkstep23(f, n, x, y, h, yh, dy);

	sum=0;
	for(int i=0; i<n; i++){
		sum += dy[i]*dy[i];}
	err = sqrt(sum);

	sum=0;
	for(int i=0; i<n; i++){
		sum += yh[i]*yh[i];}
	normy = sqrt(sum);

	tol = (normy*eps+acc)*sqrt(h/(b-a));

	if(err<tol){
		x += h;
		fprintf(list, "%20g", x);
		for(int i=0; i<n; i++){
			y[i]=yh[i];
			fprintf(list, "%20g", y[i]);
		}
		fprintf(list, "\n");
	}
	if(err>0) h*=pow(tol/err, 0.25)*0.95;
	else h*=2;
	}

	fclose(list);
}

	

