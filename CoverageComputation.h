#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

double * CoverageComputation(int *dBM,int Tbs,int Ntp,int *Xbs,int *Ybs,int *Xtp,int *Ytp){
	int i,j;
	double *coverage=malloc(Tbs*sizeof(double));

	ArrayInitialization(coverage,Tbs);

	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			if(*(dBM+i*Ntp+j)==1){
				//printf("%g %g %g %g\n",*(Xbs+i),*(Xtp+j),*(Ybs+i),*(Ytp+j));
				if(sqrt(pow((*(Xbs+i)-*(Xtp+j)),2)+pow((*(Ybs+i)-*(Ytp+j)),2))>*(coverage+i))
					*(coverage+i)=sqrt(pow((*(Xbs+i)-*(Xtp+j)),2)+pow((*(Ybs+i)-*(Ytp+j)),2));
					}
	return coverage;
}

ArrayInitialization(double *ptr, int Tbs){
       int i;
       for(i=0;i<Tbs;i++)
	*(ptr+i)=0;
}
