#include<stdio.h>
#include<math.h>
#include<omp.h>

double f(double x){
	return x*exp(2*x);
}

int main(void){
	int i = 0,n=100000000;
	double a=0, b=3, h=0, x=0, sum=0, integral=0;

/*start time*/
	double start_time=omp_get_wtime();
/*begin simpson rule*/
	#pragma omp parallel for shared(n,a,b) private(i,x) reduction(+:sum)
		for(i=1; i<n; i++){
			h=(b-a)/n;
			x=a+i*h;
			if(i%2==0){
				sum=sum+2*f(x);
			}
			else{
				sum=sum+4*f(x);
			}
		}
	double end_time=omp_get_wtime();
	printf("Time: %lf seconds \n\n",end_time-start_time);
	integral=(h/3)*(f(a)+f(b)+sum);
	printf("\n the integral is : %lf\n ",integral);
}
