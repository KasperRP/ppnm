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

	
	
	return 0;
	}

	
