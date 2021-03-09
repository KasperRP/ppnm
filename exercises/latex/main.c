#include<stdio.h>
#include<math.h>
double ex(double x);

int main(){
	int xmin=-2, xmax=2;
	double x;
	for(x=xmin;x<=xmax;x+=1.0/8)
		printf("%g %g %g\n",x,ex(x),exp(x));
return 0;
}
