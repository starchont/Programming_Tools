#include<stdio.h>
#include<omp.h>
#include<math.h>

#define N 1000

int main(void){
	int nThreads, myid, i, iStart, iEnd;
	double sum_local, sum;
	int A[N];
	sum = 0;
	for(i=0;i<N;i++){
		A[i]=0;
	}
	#pragma omp parallel shared (A,nThreads,sum)\
						private(i,iStart,iEnd,sum_local)\
						default(none)
	{
		int myid = omp_get_thread_num();
		int nThreads = omp_get_num_threads();
		iStart = myid*N/nThreads;
		iEnd = (myid+1)*N/nThreads;
		if (myid==nThreads-1){
			iEnd=N;
		}
		for(i=iStart+1;i<iEnd;i++){
			sum_local+=(3+2*i)/pow(2,i);
		}
		#pragma omp critical
		{
			sum+=sum_local;
		}
		sum_local=0;
	}
	printf("summation = %lf\n",sum);
	return 0;
}
