#include<stdio.h>
#include<limits.h>
#include<stdlib.h>
#include<assert.h>
#include<float.h>
int main() {

	int i=INT_MAX-10;
       	while(i+1>i){i++;}
	printf("INT_MAX = %i\n", INT_MAX);
	printf("\n my max int = %i\n",i); // NOTE: Code not from exercise but from examples on webpage



	for(int i=INT_MAX-10; (i+1)>i; i++){}
	printf("\n my max int = %i\n",i);

	 do{
		printf("\n my max int = %i\n",i);} while(i+1>i);

 	int j= INT_MIN + 10;
	while(j-1<j){j--;}
	printf("\nINT_MIN = %i\n",INT_MIN);
	printf("\n my min int = %i\n",j);
	
	for(int j=INT_MIN+10; (j-1)<j; j--){}
	printf("\n my min int = %i\n",j);

	do{
		printf("\n my min int = %i\n",j);} while(j-1<j);

	double xd=1;
	while(1+xd!=1){xd/=2;} xd*=2;
	printf("\nDBL_EPSILON =%g\n",DBL_EPSILON); 
	printf("\nmy dbl eps= % g\n",xd);	
	float xf=1;
	while(1+xf!=1){xf/=2;} xf*=2;
	printf("\nFLT_EPSILON = %g\n",FLT_EPSILON);
	printf("\nmy flt eps = %g\n",xf);
	long double xl=1;
	while(1+xl!=1){xl/=2;} xl*=2;
	printf("\nLDBL_EPSILON = %Lg\n",LDBL_EPSILON);
	printf("\nmy ldbl eps = %Lg\n",xl);
	
	double ed; for(ed=1; 1+ed!=1;ed/=2){} ed*=2;
	printf("\n my dbl eps = %g\n",ed);
	float ef; for(ef=1; 1+ef!=1; ef/=2){} ef*=2;
	printf("\n my flt eps = %g\n",ef);
	long double el; for(el=1; 1+el!=1; el/=2){} el*=2;
	printf("\n my lsbl eps = %Lg\n",el);

// Opgave2
	int max=INT_MAX/3;
	float sum_up_float = 0;
	float sum_down_float = 0;
	
	int n=0;
	for(n=1; n<max+1; n++){
		sum_up_float+=1.f/n;
	}
	for(n=max;n>0; n--){
		sum_down_float+=1.f/n;
	}	
	printf("sum_up_float =%f\n",sum_up_float);
	printf("sum_down_float=%f\n",sum_down_float);
	
	double sum_up_double = 0;
	double sum_down_double =0;
	for(n=1; n<max+1; n++){
		sum_up_double+=1.0/n;}
	for(n=max;n>0;n--){
		sum_down_double+=1.0/n;}
	printf("sum_up_double =%f\n",sum_up_double);
	printf("sum_down_double = %f\n", sum_down_double);

	return 0;
	}

	
