#include<math.h>
#include<assert.h>
#include<stdio.h>

// Functions to integrate in 1D

double quad(
	double f(double x),
	double a, 
	double b,
	double acc,
	double eps,
	double f2,
	double f3,
	int nrec){
assert(nrec<1000000); // max number of recurrings
double f1 = f(a+(b-a)/6);
double f4 = f(a+5*(b-a)/6);
double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // Higher order quadrature
double q = (f1+f4+f2+f3)/4*(b-a); // Lower order
double tolerance = acc+eps*fabs(Q);
double error = fabs(Q-q);

if(error < tolerance) {return Q;} // We accept if error is below tolerance
else{
	double Q1 = quad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1);
	double Q2 = quad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1);
	return Q1+Q2;} // If not accepted, we split the interval into two and recur once again
}

double integrate(
		double f(double x),
		double a, 
		double b,
		double acc,
		double eps){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	int nrec = 0;
	return quad(f, a, b, acc, eps, f2, f3, nrec);
}

// 2D now

double quad2D(
	double F(double x, double y),
	double a,
	double b,
	double d(double x),
	double u(double x),
	double acc,
	double eps,
	double F2, 
	double F3,
	int nrec){
assert(nrec<1000000);
double X1 = a+(b-a)/6;
double X4 = a+5*(b-a)/6;
double f1(double y){ return F(X1, y);}
double f4(double y){ return F(X4, y);}
double F1 = integrate(f1, d(X1), u(X1), acc, eps);
double F4 = integrate(f4, d(X4), u(X4), acc, eps);
double I = (2*F1+F2+F3+2*F4)/6*(b-a);
double i = (F1+F4+F2+F3)/4*(b-a);
double tolerance = acc+eps*fabs(I);
double error = fabs(I-i);
if(error < tolerance) { return I;}
else{
	double I1 = quad2D(F, a, (a+b)/2, d, u, acc/sqrt(2), eps, F1, F2, nrec+1);
	double I2 = quad2D(F, (a+b)/2, b, d, u, acc/sqrt(2), eps, F3, F4, nrec+1);
	return I1+I2;}
}

double integrate2D(
		double F(double x, double y),
		double a,
		double b,
		double d(double x),
		double u(double x),
		double acc,
		double eps){
double X2 = a+2*(b-a)/6;
double X3 = a+4*(b-a)/6;
double f2(double y) {return F(X2,y);}
double f3(double y) {return F(X3,y);}
double F2 = integrate(f2, d(X2), u(X2), acc, eps);
double F3 = integrate(f3, d(X3), u(X3), acc, eps);
int nrec = 0;
return quad2D(F, a, b, d, u, acc, eps, F2, F3, nrec);
}
