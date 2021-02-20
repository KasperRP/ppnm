#include<stdio.h>
#include<math.h>

int main() {
	double x; 
	int items;

	FILE* inputfile=fopen("input.txt","r");
	FILE* outputfile=fopen("out.file.txt","w");

	do{
		items=fscanf(inputfile,"%lg",&x);
		fprintf(outputfile,"x=%g sin(x)=%g cos(x)=%g\n",x,sin(x),cos(x));
	} while(items != EOF);

	fclose(inputfile);
	fclose(outputfile);
return 0;
}
