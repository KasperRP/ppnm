#include<math.h>
#include<assert.h>
#include<stdio.h>

double adaptint(
	double f(double x), 
	double a, 
	double b, 
	double acc, 
	double eps, 
	double f2, 
	double f3, 
	int nrec){
assert(nrec<10000000); // max number of recurrings
double f1 = f(a+(b-a)/6);
double f4=f(a+5*(b-a)/6);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a); // Higher order quadrature
double q=(f1+f4+f2+f3)/4*(b-a); // Lower order
double tolerance = acc+eps*fabs(Q);
double error = fabs(Q-q);

if(error < tolerance){ return Q;} // We accept if error is below tolerance
else{
	double Q1=adaptint(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1);
	double Q2=adaptint(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1);
	return Q1+Q2;} // If not accepted, we split the interval into two and recur once again
}

double integrate(
		double f(double x),
		double a,
		double b, 
		double acc, 
		double eps){
	double f2=f(a+2*(b-a)/6);
	double f3=f(a+4*(b-a)/6);
	int nrec=0;
	return adaptint(f, a, b, acc, eps, f2, f3, nrec);
}


double CCintegrate( // Clenshaw-Curtis integral
		double f(double),
		double a,
		double b,
		double acc,
		double eps){
	
	double F(double t){
		return f((b-a)*cos(t)/2+(b+a)/2)*sin(t)*(b-a)/2;}

	

	return integrate(F,0,M_PI,acc,eps);
}




//static double A;
//static double B;
//double F(double f(double), double t){
//	return f((A+B)/2+(A-B)/2*cos(t))*sin(t)*(B-A)/2;
//}

//double CCadaptint(
//		double f(double),
//		double a,
//		double b,
//		double acc,
//		double eps, 
//		double f2,
//		double f3,
//		int nrec){
//	assert(nrec<1000000);
//	double f1=F(f,a+(b-a)/6);
//	double f4=F(f,5*(b-a)/6);
//	double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
//	double q=(f1+f4+f2+f3)/4*(b-a);
//	double tol=acc+eps*fabs(Q);
//	double err=fabs(Q-q);
//	if(err<tol) {return Q;}
//	else{
//		double Q1=CCadaptint(f,a,(a+b)/2,acc/sqrt(2),eps,f1,f2,nrec+1);
//		double Q2=CCadaptint(f,(a+b)/2,b,acc/sqrt(2),eps,f3,f4,nrec+1);
//		return Q1+Q2;}
//}

//double CCintegrate(
//		double f(double),
//		double a,
//		double b,
//		double acc,
//		double eps){
//	A=a; B=b;
//	a=0; b=M_PI;
//	double f2=F(f,a+2*(b-a)/6);
//	double f3=F(f,a+4*(b-a)/6);
//	int nrec=0;
//	return CCadaptint(f,a,b,acc,eps,f2,f3,nrec);
//}


