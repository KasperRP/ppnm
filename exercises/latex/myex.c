#include<math.h>

double ex(double x){
	if(x<0) return 1/ex(-x);
	if(x>1./8) return pow(ex(x/2),2);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}
