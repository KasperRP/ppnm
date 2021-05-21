#include<math.h>
#include<assert.h>
#include<stdio.h>



double quad(
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
	double Q1=quad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1);
	double Q2=quad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1);
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
	return quad(f, a, b, acc, eps, f2, f3, nrec);
}


// Functions used for variable transformations

static double A;
static double B; // Not accessible from other files

double F(double f(double), double t){ // Transformed function for Clenshaw-Curtis

	return f((A+B)/2+(A-B)/2*cos(t))*sin(t)*(B-A)/2;
}

// Transformed functions for infinite integral

double Fint(double f(double), double t){ // For both limits infinte
	
	return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t));}

double Gint(double f(double), double t){

	return f(A+t/(1-t))/((1-t)*(1-t));}

double Hint(double f(double), double t){

	return f(B+t/(1+t))/((1+t)*(1+t));}

// Clenshaw-Curtis implementation:

double CCquad( // Just like quad but use F instead of f
	double f(double x),
	double a,
	double b,
	double acc,
	double eps,
	double f2,
	double f3,
	int nrec){
	assert(nrec<1000000); // max number of recurrings
	double f1 = F(f, a+(b-a)/6);
	double f4 = F(f, a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // Higher order quadrature
	double q = (f1+f4+f2+f3)/4*(b-a); // Lower order
	double tolerance = acc+eps*fabs(Q);
	double error = fabs(Q-q);

	if(error < tolerance) return Q;
	else{
		double Q1 = CCquad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1);
		double Q2 = CCquad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1);
		return Q1+Q2;
	}
}


double CCintegrate( // Clenshaw-Curtis integral
		double f(double),
		double a,
		double b,
		double acc,
		double eps){
	A=a; B=b;
	a=0; b=M_PI;
	double f2=F(f,a+2*(b-a)/6);
	double f3=F(f,a+4*(b-a)/6);
	int nrec=0;

	return CCquad(f,a,b,acc,eps,f2,f3,nrec);
}


//double CCintegrate(  OBS: Nested functions do not work for Linux subspace in Windows!
//		double f(double),
//		double a,
//		double b,
//		double acc,
//		double eps){

//	double F(double t){
//		return f( (a+b)/2+(a-b)/2*cos(t) )*sin(t)*(b-a)/2;}
	
//	return integrate(F, 0, M_PI, acc, eps);
//}


/*
double doubleinfquad(
		double f(double x),
		double a,
		double b,
		double acc,
		double eps,
		double f2,
		double f3,
		int nrec,
		double* error){

	assert(nrec<1000000);
	double f1 = Fint(f, a+(b-a)/6);
	double f4 = Fint(f, a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f4+f2+f3)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	double err = fabs(Q-q);
	
	if(err < tol) { error* = err; return Q;}
	else{
		double Q1 = doubleinfquad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1, error);
		double Q2 = doubleinfquad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1, error);
		return Q1+Q2;
	}
}

double upperinfquad(
		double f(double x),
		double a,
		double b,
		double acc,
		double eps,
		double f2,
		double f3,
		int nrec,
		double* error){

	assert(nrec<1000000);
	double f1 = Gint(f, a+(b-a)/6);
	double f4 = Gint(f, a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f4+f2+f3)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	double err = fabs(Q-q);

	if(err < tol) { error* = err; return Q;}
	else{
		double Q1 = upperinfquad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1, error);
		double Q2 = upperinf(quad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1, error);
		return Q1+Q2;
		}
	}

double lowerinfquad(
		double f(double x),
		double a,
		double b,
		double acc,
		double eps,
		double f2,
		double f3,
		int nrec,
		double* error){

	assert(nrec<10000000);
	double f1 = Hint(f, a+(b-a)/6);
	double f4 = Hint(f, a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f4+f2+f3)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	double err = fabs(Q-q);

	if(err < tol) {error* = err; return Q;}
	else{
		double Q1 = lowerinfquad(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1, error);
		double Q2 = lowerinfquad(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1, error);
		return Q1+Q2;
	}
}

double infintegrate(
		double f(double),
		double a,
		double b,
		double acc,
		double eps){

	if(isinf(a)==-1 && isinf(b)==1){
		
		double xmin = -1;
		double xmax = 1;
		double f2 = Fint(f, xmin+2*(xmax-xmin)/6);
		double f3 = Fint(f, xmin+4*(xmax-xmin)/6);
		int nrec = 0;
		return doubleinfquad(f, xmin, xmax, acc, eps, f2, f3, nrec, error);
	}

	if(isinf(a)==0 && isinf(b)==1){

		A=a;
		double xmin = 0;
		double xmax = 1;
		double f2 = Gint(f, xmin+2*(xmax-xmin)/6);
		double f3 = Gint(f, xmin+4*(xmax-xmin)/6);
		int nrec = 0;
		return upperinfquad(f, xmin, xmax, acc, eps, f2, f3, nrec, error);
	}

	if(isinf(a)==-1 && isinf(b)==0){

		B=b;
		double xmin = -1; // See notes for transform integrals
		double xmax = 0; 
		double f2 = Hint(f, xmin+2*(xmax-xmin)/6);
		double f3 = Hint(f, xmin+4*(xmax-xmin)/6);
		int nrec = 0;
		return lowerinfquad(f, xmin, xmax, acc, eps, f2, f3, nrec, error);
	}
}
*/


