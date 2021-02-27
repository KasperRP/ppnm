#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f(double x, void* params) {
	double z = *(double*) params;
	double f = log(x*z)/sqrt(x);
	return f;
}

int main() {
	gsl_function F;
	F.function = &f;
	double z = 1.0;
	F.params = (void*)&z;
	int limit = 999;
	gsl_integration_workspace *w;
	w = gsl_integration_workspace_alloc (limit);
	double a = 0, b = 1, acc = 1e-6, eps = 1e-6, result, error;
	gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
	gsl_integration_workspace_free(w);
	printf("result = %g\n", result);
return 0;
}
