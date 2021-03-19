#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include"binsearch.h"
int binsearch(gsl_vector* x, double z){
	int n = (*x).size;
	assert(gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid;
		else j=mid;
	}
return i; // returns the left side of interval
}
