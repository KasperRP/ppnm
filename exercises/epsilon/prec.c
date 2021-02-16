#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int equal(double a,double b,double tau,double eps){
		if(fabs(a-b)<tau || fabs(a-b)/(fabs(a)+fabs(b))<eps/2)
		{return 1;}
 		 else{return 0;}
}
	
int main(){
	printf("a = 1, b = 2, tau = 1, eps = 2 leads to %i\n",equal(1,2,1,2));
	return 0;
}
