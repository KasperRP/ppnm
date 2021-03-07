#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<pthread.h>
// see harm.c inside Lecture-8 folder

struct params {int N, count; unsigned int seed;};   // seed is for the rand_r command

void* myrand(void* arg) {               // this type of function is needed for pthread_create
	struct params* p = (struct params*) arg;
	int N = (*p).N;
	unsigned int seed = (*p).seed;
	int n=0;
	
	for(int i=0; i<=N;i++){
		double x = (double) rand_r(&seed)/RAND_MAX;
		double y = (double) rand_r(&seed)/RAND_MAX;
		if(x*x+y*y<=1) n++;}
	(*p).count = n;
	return NULL;
}

int main() {
	int N =(int) 1e8;
	pthread_t t1,t2,t3;
	struct params p1 = {.N=N/3, .count=0, .seed=1};
	struct params p2 = {.N=N/3, .count=0, .seed=2};
	struct params p3 = {.N=N/3, .count=0, .seed=3};
	pthread_create(&t1,NULL,myrand,(void*)&p1);
	pthread_create(&t2,NULL,myrand,(void*)&p2);
	pthread_create(&t3,NULL,myrand,(void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	double count = p1.count+p2.count+p3.count;
	double pi = 4*count/N;
	printf("Estimate of pi = %g\n",pi);
return 0;
}	

