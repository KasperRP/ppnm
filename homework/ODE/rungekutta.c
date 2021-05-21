#include<stdio.h>
#include<math.h>
#include<stdlib.h>

// I implement the stepper based on the Bogacki-Shampine method (see notes). 

void rkstep23(
	void f(int n, double x, double* y, double* dydx), //f(x,y)=dydx. Evaluation of f stored in dydx 
	int n, // size of vector y
	double x,
	double* yx, // y(x)
	double h, 
	double* yh, // y(x+h)
	double* dy) { // error

	double k1[n], k2[n], k3[n], k4[n];
	double yt[n];

	f(n, x, yx, k1); // Evaluates to determine k1
	for(int i=0; i<n; i++) {yt[i]=yx[i]+1./2*k1[i]*h;} 

	f(n, x+1./2*h, yt, k2);
	for(int i=0; i<n; i++) {yt[i]=yx[i]+3./4*k2[i]*h;}

	f(n, x+3./4*h, yt, k3); 
	for(int i=0; i<n; i++) {yh[i]=yx[i]+(2./9*k1[i]+1./3*k2[i]+4./9*k3[i])*h;}

	f(n, x+h, yh, k4);
	for(int i=0; i<n; i++){
		yt[i]=yx[i]+(7./24*k1[i]+1./4*k2[i]+1./3*k3[i]+1./8*k4[i])*h;
		dy[i]=yh[i]-yt[i];
	}
}


// Next I create the driver

void driver(
		void f(int n, double x, double* y, double* dydx),
		int n,
		double a,
		double* ya,
		double b,
		//double* y, // y starts as y(a) and ends as y(b)
		double* yb,
		double h, 
		double acc,
		double eps,
		char* path // Storing the path (B)
	   )
{

	FILE* list = fopen(path, "w");
	double x;
	double y[n];
       	double sum;
        double	err;
        double normy;
        double tol;
        double yh[n];
        double	dy[n];

	// First step
	x=a;
	fprintf(list, "%20g", x);
	for(int i=0; i<n; i++){
		y[i]=ya[i];
		fprintf(list, "%20g", y[i]);
	}
	fprintf(list, "\n");

	// starting while-loop
	while(x<b){

		if(x+h>b) h=b-x;
	
	// make step
	rkstep23(f, n, x, y, h, yh, dy);
	
	// Calculating total error err from local error dy^2
	
	sum=0;
	for(int i=0; i<n; i++){
		sum += dy[i]*dy[i];}
	err = sqrt(sum);

	// Calculating norm of estimated y(x+h)
	
	sum=0;
	for(int i=0; i<n; i++){
		sum += yh[i]*yh[i];}
	normy = sqrt(sum);

	// Estimate tolerance
	tol = (normy*eps+acc)*sqrt(h/(b-a));

	if(err<tol){ // if err<tau we accept
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
	}  // end of while loop
	
	for(int i=0; i<n; i++){
		yb[i]=yh[i];
	}
	

	fclose(list);
}




