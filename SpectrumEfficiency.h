#include<stdio.h> 
#include<time.h>
#include<stdlib.h>
#include"S2Dquicksort4.h"
#include"Squicksort2.h"

#define SBP 8.00	//Power for BS in sleep mode
#define ABP 68.73	//Basic power for BS in active mode
//extern int MNbs;
/*extern int Ntp;
extern int MP;
extern int DSt;
extern double BP;
extern double P[][];
extern int DS[];*/

double * SCoverageFinder(int *dBM, double *D, int Tbs, int Ntp){
	int i,j;
	double *coverage=malloc(Tbs*sizeof(double));
	SArrayInitialization(coverage,Tbs);
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			if(*(dBM+i*Ntp+j)==1)
				if(((*(D+(i+1)*(Ntp+1)+j+1))/20)>coverage[i])
					coverage[i]=(*(D+(i+1)*(Ntp+1)+j+1))/20;
					//printf("%g\n",coverage[i]);
	return coverage;
}

SdBMprinter(int *dBM, int Ntp, int Tbs, int *DS, double *P){
               int i,j;
               for(i=0;i<Tbs;i++)
               printf("       BS%3d",i+1);
               printf("\n");
               
               for(i=0;i<Ntp;i++){
                                  printf("MS%4d%4d",i+1,*(dBM+i));
                                  for(j=1;j<Tbs;j++)
                                  printf("      %6d",*(dBM+j*Ntp+i));
                                  printf("\n");
                                  
                                  printf(" Power%10.5f",*(P+(Ntp+1)+i+1));
                                  for(j=1;j<Tbs;j++)
                                  printf("  %10.5f",*(P+(j+1)*(Ntp+1)+i+1));
                                  printf("\n");
                                  
                                  printf(" DS %6d",*(DS+i+1));
                                  for(j=1;j<Tbs;j++)
                                  printf("      %6d",*(DS+i+1));
                                  printf("\n");
                                  }
                                  }

Stemp_dBM_Initiator(int *ptr, int Ntp){	//initiate the array of temp_dBM for the new heuristic algorithm
	int i;
	for(i=0;i<Ntp;i++)
		*(ptr+i)=-1;
}

SArrayInitialization(double *ptr, int Tbs){
       int i;         
       for(i=0;i<Tbs;i++)
                          *(ptr+i)=0;

       
}
SArrayInitialization1(int *ptr, int Tbs){
       int i;         
       for(i=0;i<Tbs;i++)
                          *(ptr+i)=0;

       
}     
SArrayInitialization2(int *ptr, int Tbs, int Ntp){
       int i,j;         
       for(i=0;i<Tbs;i++)
       for(j=0;j<Ntp;j++)
       *(ptr+i*Ntp+j)=0;

       
}           
SArrayInitialization3(double *ptr, int Tbs, int Ntp){
	int i,j;
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			*(ptr+i*Ntp+j)=0.0;
}
Sheuristic(int Tbs, double *P, int Ntp, int MP, int DSt, double BP, int *DS, int *Xbs, int *Ybs, int *Xtp, int *Ytp, double *ptrD, double *BR){
                
       //double *SortedPowerRequirement;  
       //int *SortedPowerRequirementIndex;
	double UnsortedObjective[Tbs][Ntp];
	SArrayInitialization3(UnsortedObjective[0],Tbs,Ntp);
	double SortedObjective[Tbs][Ntp];
	int *SortedObjectiveIndex;

	double UnsortedCumulativeObjective[Tbs];
	SArrayInitialization(UnsortedCumulativeObjective,Tbs);
	double SortedCumulativeObjective[Tbs];
	double ObjectiveValue=0.0;
	int *SortedCumulativeObjectiveIndex;
       
       double UnsortedCumulativePower[Tbs];
       SArrayInitialization(UnsortedCumulativePower,Tbs);
       //double SortedCumulativePower[Tbs];
       double TotalPower=0.0;
       //int *SortedCumulativePowerIndex;
       
       int CumulativeDS[Tbs];
       SArrayInitialization1(CumulativeDS,Tbs);
       int TotalDSs=0;

	double CumulativeBR[Tbs];
	SArrayInitialization(CumulativeBR,Tbs);
	double TotalBR=0.0;
	double NOV=0.0;

	int BStatus[Tbs];
	SArrayInitialization1(BStatus,Tbs);
                            
       int dBM[Tbs][Ntp];
       SArrayInitialization2(dBM[0],Tbs,Ntp);
       int CumulativeServingMSs[Tbs];
       SArrayInitialization1(CumulativeServingMSs,Tbs);
       int TotalServedMSs=0;
       
       int i,j,k,l;
       //double Max=-10000000000.0;
       //int IndexMax=-1;

       double *coverage;
           
       clock_t start_tick, end_tick;
       double elapsed;
       
       start_tick=clock();

	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++){
			UnsortedObjective[i][j]=(*(BR+(j+1)))/(*(DS+(j+1)))/(*(P+(i+1)*(Ntp+1)+j+1));
			SortedObjective[i][j]=(*(BR+(j+1)))/(*(DS+(j+1)))/(*(P+(i+1)*(Ntp+1)+j+1));
			//printf("BS%d MS%d: %lf \n",i+1,j+1,UnsortedObjective[i][j]);
 			//printf("BR:%lf DS:%d Power:%lf\n",*(BR+(j+1)),*(DS+(j+1)),*(P+(i+1)*(Ntp+1)+j+1));
		}
	//SortedObjective=Unsorted_to_Sorted(UnsortedObjective,Tbs,Ntp);
	SortedObjectiveIndex=Sq_sort_2D(*SortedObjective,0,Tbs-1,Ntp,Tbs);
       //SortedPowerRequirement=SortedPowerRequirementAssignment(P,Tbs,Ntp);
       		/*	
                         for(i=0;i<Tbs;i++)
                         for(j=0;j<Ntp;j++)
                         printf("BS%d MS%d: %lf \n",i+1,j+1,SortedObjective[i][j]);
                         */
       //SortedPowerRequirementIndex=q_sort_2D(SortedPowerRequirement,0,Tbs-1,Ntp,Tbs);
       
      
       /*for(k=0;k<Ntp;k++){
       for(j=0;j<Tbs;j++)
       printf("Sorted Objective from BS %d to MS %d: %lf \n",*(SortedObjectiveIndex+j*Ntp+k)+1,k+1,SortedObjective[j][k]);
       printf("\n");
       }*/
       
       //spectrum algorithm goes from here
       int temp_dBM[Ntp];
       Stemp_dBM_Initiator(temp_dBM,Ntp);
	double Min=1000.0;
	int IndexMin=-1;	

       for(i=0;i<Ntp;i++){	//each MS is assigned to the 1st BS choice
	*(temp_dBM+i)=*(SortedObjectiveIndex+i+Ntp*(Tbs-1));
	CumulativeServingMSs[*(SortedObjectiveIndex+i+Ntp*(Tbs-1))]++;
	TotalServedMSs++;
	*(UnsortedCumulativePower+*(SortedObjectiveIndex+i+Ntp*(Tbs-1)))+=*(P+(*(SortedObjectiveIndex+i+Ntp*(Tbs-1))+1)*(Ntp+1)+i+1);
	*(CumulativeDS+*(SortedObjectiveIndex+i+Ntp*(Tbs-1)))+=*(DS+i+1);
		*(CumulativeBR+*(SortedObjectiveIndex+i+Ntp*(Tbs-1)))+=*(BR+i+1);
		//printf("%lf %lf\n",*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-1))),SortedObjective[Tbs-1][i]);
		*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-1)))+=SortedObjective[Tbs-1][i];
		//printf("%lf\n",*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-1))));
	BStatus[*(SortedObjectiveIndex+i+Ntp*(Tbs-1))]=1;
}
int temp=0;
for(i=0;i<Tbs;i++){
	printf("%d\n", BStatus[i]);
	 //printf("Unsorted Cumulative Objective of BS %d: %lf \n",i+1,*(UnsortedCumulativeObjective+i));
}
/*if(TotalServedMSs!=temp){
printf("0. Inconsistent numbers! %d %d %d %d %d %d",TotalServedMSs,CumulativeServingMSs[0],CumulativeServingMSs[1],CumulativeServingMSs[2],CumulativeServingMSs[3],temp);
exit(1);}*/

for(j=0;j<Tbs;j++){
	while((*(UnsortedCumulativePower+j)>MP)||(*(CumulativeDS+j)>DSt)){
		for(i=0;i<Ntp;i++){
			if(*(temp_dBM+i)==j&&UnsortedObjective[*(temp_dBM+i)][i] < Min){
				Min=UnsortedObjective[*(temp_dBM+i)][i];
				IndexMin=i;
			}
		}
		if(IndexMin!=-1){
			*(UnsortedCumulativePower+*(temp_dBM+IndexMin))-=*(P+(*(temp_dBM+IndexMin)+1)*(Ntp+1)+IndexMin+1);
			*(UnsortedCumulativeObjective+*(temp_dBM+IndexMin))-=UnsortedObjective[*(temp_dBM+IndexMin)][IndexMin];
			*(CumulativeDS+*(temp_dBM+IndexMin))-=*(DS+IndexMin+1);
			*(CumulativeBR+*(temp_dBM+IndexMin))-=*(BR+IndexMin+1);
			CumulativeServingMSs[*(temp_dBM+IndexMin)]--;
			TotalServedMSs--;
			*(temp_dBM+IndexMin)=-1;
		}
	BStatus[j]=2;
	IndexMin=-1;
	Min=1000.0;
	}	
}

if(TotalServedMSs<ceil(Ntp*(1-BP))){
	for(i=0;i<Ntp;i++){
		if(*(temp_dBM+i)==-1){
			*(UnsortedCumulativePower+*(SortedObjectiveIndex+i+Ntp*(Tbs-2)))+=*(P+(*(SortedObjectiveIndex+i+Ntp*(Tbs-2))+1)*(Ntp+1)+i+1);
			*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+i+Ntp*(Tbs-2)))+=SortedObjective[Tbs-2][i];
			*(CumulativeDS+*(SortedObjectiveIndex+i+Ntp*(Tbs-2)))+=*(DS+i+1);
			*(CumulativeBR+*(SortedObjectiveIndex+i+Ntp*(Tbs-2)))+=*(BR+i+1);
			*(temp_dBM+i)=*(SortedObjectiveIndex+i+Ntp*(Tbs-2));
			CumulativeServingMSs[*(SortedObjectiveIndex+i+Ntp*(Tbs-2))]++;
			TotalServedMSs++;
			BStatus[*(SortedObjectiveIndex+i+Ntp*(Tbs-2))]=1;
		}
	}
}

for(j=0;j<Tbs;j++){
	while((*(UnsortedCumulativePower+j)>MP)||(*(CumulativeDS+j)>DSt)){
		for(i=0;i<Ntp;i++){
			if((*(temp_dBM+i)==j)&&UnsortedObjective[*(temp_dBM+i)][i]<Min){
				Min=UnsortedObjective[*(temp_dBM+i)][i];
				IndexMin=i;
			}
		}
		if(IndexMin!=-1){
			*(UnsortedCumulativePower+*(temp_dBM+IndexMin))-=*(P+(*(temp_dBM+IndexMin)+1)*(Ntp+1)+IndexMin+1);
			*(UnsortedCumulativeObjective+*(temp_dBM+IndexMin))-=UnsortedObjective[*(temp_dBM+IndexMin)][IndexMin];
			*(CumulativeDS+*(temp_dBM+IndexMin))-=*(DS+IndexMin+1);
			*(CumulativeBR+*(temp_dBM+IndexMin))-=*(BR+IndexMin+1);
			CumulativeServingMSs[*(temp_dBM+IndexMin)]--;
			TotalServedMSs--;
			*(temp_dBM+IndexMin)=-1;
		}
	BStatus[j]=2;
	IndexMin=-1;
	Min=1000.0;
	}	
}
printf("TotalServedMSs %d\n",TotalServedMSs);
for(i=0;i<Tbs;i++)
	printf("CumulativeServingMSs %d\n",CumulativeServingMSs[i]);
for(i=0;i<Tbs;i++){
		if(*(UnsortedCumulativePower+i)>0){
			*(UnsortedCumulativePower+i)+=ABP;
			//*(UnsortedCumulativeObjective+i)-=ABP;
		}
		else{
			*(UnsortedCumulativePower+i)+=SBP;
			//*(UnsortedCumulativeObjective+i)-=SBP;
		}
}
for(i=0;i<Tbs;i++)
		printf("BS%d: %lf \n",i+1,*(UnsortedCumulativeObjective+i));
for(i=0;i<Tbs;i++)
       *(SortedCumulativeObjective+i)=*(UnsortedCumulativeObjective+i);
       
       printf("\n");
       SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);/*
       for(i=0;i<Tbs;i++)
       printf("Unsorted Cumulative Objective of BS %d: %lf \n",i+1,*(UnsortedCumulativeObjective+i));
       printf("\n");
       for(i=0;i<Tbs;i++)
       printf("%d. Sorted Cumulative Objective of BS %d: %lf \n",i+1,*(SortedCumulativeObjectiveIndex+i)+1,*(SortedCumulativeObjective+i));
       printf("\n");
       for(i=0;i<Tbs;i++)
       printf("%d ",*(SortedCumulativeObjectiveIndex+i)+1);
       printf("\n");

	for(i=0;i<Tbs;i++)
       printf("Unsorted Cumulative Power of BS %d: %lf \n",i+1,*(UnsortedCumulativePower+i));
       printf("\n");
       */
       
	double difference;
	int RemovedUsers=0;
  	//for(i=Tbs-1;i>-1;i--){
	for(i=0;i<Tbs;i++){
		//while(difference>0&&difference<ABP-SBP)
			if(BStatus[*(SortedCumulativeObjectiveIndex+i)]==1){
				difference=0.0;
				RemovedUsers=0;
				int *AdditionDS=malloc(Tbs*sizeof(int));
				double *AdditionPower=malloc(Tbs*sizeof(double));
				SArrayInitialization(AdditionPower,Tbs);
				SArrayInitialization1(AdditionDS,Tbs);
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i))
						for(k=Tbs-1;k>-1;k--)
							if(BStatus[*(SortedObjectiveIndex+k*Ntp+j)]==1&&*(SortedObjectiveIndex+k*Ntp+j)!=*(SortedCumulativeObjectiveIndex+i)&&*(AdditionDS+*(SortedObjectiveIndex+k*Ntp+j))+*(DS+j+1)+(*(CumulativeDS+*(SortedObjectiveIndex+k*Ntp+j)))<=DSt&&*(AdditionPower+*(SortedObjectiveIndex+k*Ntp+j))+*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)+(*(UnsortedCumulativePower+*(SortedObjectiveIndex+k*Ntp+j)))-ABP<=MP){
								difference+=*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)-*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
								
								RemovedUsers++;
								//printf("Removalbe User Number: %d\n",RemovedUsers);
								*(AdditionDS+*(SortedObjectiveIndex+k*Ntp+j))+=*(DS+j+1);
								*(AdditionPower+*(SortedObjectiveIndex+k*Ntp+j))+=*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1);
								k=-1;
							}
				}
				free(AdditionPower);
				free(AdditionDS);
				if(difference>0&&difference<ABP-SBP&&RemovedUsers==*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))){
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i))
							for(k=Tbs-1;k>-1;k--)
								if(BStatus[*(SortedObjectiveIndex+k*Ntp+j)]==1&&*(SortedObjectiveIndex+k*Ntp+j)!=*(SortedCumulativeObjectiveIndex+i)&&(*(CumulativeDS+*(SortedObjectiveIndex+k*Ntp+j))+*(DS+j+1))<=DSt&&(*(UnsortedCumulativePower+*(SortedObjectiveIndex+k*Ntp+j))+*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)-ABP)<=MP){
									//if((*(CumulativeDS+*(SortedObjectiveIndex+k*Ntp+j))+*(DS+j+1))<=DSt&&(*(UnsortedCumulativePower+*(SortedObjectiveIndex+k*Ntp+j))+*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1)-ABP)<=MP){
										*(UnsortedCumulativePower+*(SortedObjectiveIndex+k*Ntp+j))+=*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1);
										*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
										*(CumulativeDS+*(SortedObjectiveIndex+k*Ntp+j))+=*(DS+j+1);
										*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
										if((*(UnsortedCumulativePower+*(SortedObjectiveIndex+k*Ntp+j))-ABP)==MP||(*(CumulativeDS+*(SortedObjectiveIndex+k*Ntp+j))==DSt))
											BStatus[*(SortedObjectiveIndex+k*Ntp+j)]=2;
										*(UnsortedCumulativeObjective+*(SortedObjectiveIndex+k*Ntp+j))+=(*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedObjectiveIndex+k*Ntp+j)+1)*(Ntp+1)+j+1));
										*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=(*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1));
										*(CumulativeBR+*(SortedObjectiveIndex+k*Ntp+j))+=*(BR+j+1);
										*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);
										*(CumulativeServingMSs+*(SortedObjectiveIndex+k*Ntp+j))+=1;
										*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
										//printf("BS%d has users %d\n",*(SortedCumulativeObjectiveIndex+i)+1,*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i)));
										*(temp_dBM+j)=*(SortedObjectiveIndex+k*Ntp+j);
										k=-1;
									//}
									
								}
					if(*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))==0){
						BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
						*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+=SBP-ABP;
						//*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))+=ABP-SBP;
					}
					else{
						printf("1. Program Stops! BS%d serves %d users instead of 0.\n",*(SortedCumulativeObjectiveIndex+i)+1,*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i)));
						exit(1);
					}
					for(l=0;l<Tbs;l++)
						if(*(UnsortedCumulativePower+l)-ABP>MP||*(CumulativeDS+l)>DSt){
							printf("2. Program Stops!\n");
							exit(1);
							//BStatus[l]=2;
						}
					for(l=0;l<Tbs;l++)
						*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
					SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
					i=0;
				}
			}
	}
	for(i=0;i<Tbs;i++)
		printf("%d ",BStatus[i]);
	printf("\n");
	for(l=0;l<Tbs;l++)
		printf("BS:%d Obj:%lf\n",*(SortedCumulativeObjectiveIndex+l)+1,*(SortedCumulativeObjective+l));
	printf("\n");


double NewTransmissionPower=0.0;

double NewSE;
double MaxIncrease;
int NewBS=-1;
int UserCutOff;
double MaxP=0.0;
int IndexMaxP=-1;
for(i=0;i<Tbs;i++){
	NewBS=-1;
	MaxIncrease=0.0;
	UserCutOff=-1;
// 	NewPower=0.0;
	//printf("%d. BS status:%d\n",i,BStatus[*(SortedCumulativeObjectiveIndex+i)]);
	int b,u;
	double se=0.0;
	double p=0.0;
	double obj=0.0;
	for(u=0;u<Ntp;u++)
		if(*(temp_dBM+u)!=-1)
			se+=(*(BR+u+1))/(*(DS+u+1));
	for(b=0;b<Tbs;b++)
		p+=*(UnsortedCumulativePower+b);
	obj=se/p;
	if(BStatus[*(SortedCumulativeObjectiveIndex+i)]==1||BStatus[*(SortedCumulativeObjectiveIndex+i)]==2){
// 		printf("BS status:%d\n",BStatus[*(SortedCumulativeObjectiveIndex+i)]);
		for(k=Tbs-1;k>-1;k--){
			if(*(SortedCumulativeObjectiveIndex+k)!=*(SortedCumulativeObjectiveIndex+i)){
				NewTransmissionPower=0.0;
// 				NewObjective=0.0;
				NewSE=0.0;
				//printf("Old BS:%d, New BS:%d\n",*(SortedCumulativePowerIndex+i)+1,*(SortedCumulativePowerIndex+k)+1);
				if((BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)||(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0)){
// 					if(((*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k)))+(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))))<=DSt)
						for(j=0;j<Ntp;j++){
							if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
								NewTransmissionPower+=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
// 								NewSE+=(*(BR+j+1))/(*(DS+j+1));
// 								NewObjective+=(*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1));
								if((*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1))!=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j]){
									printf("NewObjective:%lf UnsotedObjecitve:%lf\n",(*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1)),UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j]);
									exit(1);
								}
							}
						}
				}
 				printf("To New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",*(SortedCumulativeObjectiveIndex+k)+1,NewTransmissionPower,*(SortedCumulativeObjectiveIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i)));		
				if(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0){
					/*for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i))*/
					if(NewTransmissionPower<=MP){
						NewSE=se/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))+NewTransmissionPower+ABP+SBP);
						if(NewSE-obj>MaxIncrease){
							UserCutOff=0;
							NewBS=*(SortedCumulativeObjectiveIndex+k);
							MaxIncrease=NewSE-obj;
							printf("Users in BS%d is gonna be reassigned to BS%d with an additional power level of %lf\n",*(SortedCumulativeObjectiveIndex+i)+1,NewBS+1,NewTransmissionPower);
						}
					}
					else if(NewTransmissionPower>MP){
						int MergedSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
						double *MergedPower=calloc(MergedSize,sizeof(double));
 						double *MergedSE=calloc(MergedSize,sizeof(double));
						int *MergedIndex;
						int ExceedUser=TotalServedMSs-ceil(Ntp*(1-BP));
						double ExceedPower=0.0;
 						double ExceedSE=0.0;
						int m=0;
						if(ExceedUser>0){
							for(j=0;j<Ntp;j++)
								if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
									*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
									*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
									m++;
								}
							for(j=0;j<Ntp;j++)
								if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
									*(MergedPower+m)=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
									*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
									m++;
								}
							MergedIndex=Sq_sort(MergedPower,0,MergedSize-1,MergedSize);
							j=MergedSize-1;
							if(MergedSize>ExceedUser)
								while(j>MergedSize-ExceedUser-1){
									ExceedPower+=MergedPower[j];
									ExceedSE+=MergedSE[*(MergedIndex+j)];
									NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))+ABP+SBP+NewTransmissionPower-ExceedPower);
									if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-SBP-ExceedPower<=MP){
										NewBS=*(SortedCumulativeObjectiveIndex+k);
										MaxIncrease=NewSE-obj;
										UserCutOff=1;
										j=-1;
									}
									else
										j--;
								}
							else if(MergedSize<=ExceedUser)
								while(j>-1){
									ExceedPower+=MergedPower[j];
									ExceedSE+=MergedSE[*(MergedIndex+j)];
									NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))+ABP+SBP+NewTransmissionPower-ExceedPower);
									if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-SBP-ExceedPower<=MP){
										NewBS=*(SortedCumulativeObjectiveIndex+k);
										MaxIncrease=NewSE-obj;
										UserCutOff=1;
										j=-1;
									}
									else
										j--;
								}
						}
						free(MergedPower);
						free(MergedSE);
// 						delete MergedIndex;
					}
				}
				else if(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1){
					if(NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))<=MP+ABP){
						if(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))<=DSt){
							NewSE=se/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+NewTransmissionPower+SBP);
							if(NewSE-obj>MaxIncrease){
								UserCutOff=0;
								NewBS=*(SortedCumulativeObjectiveIndex+k);
								MaxIncrease=NewSE-obj;
							}
						}
						else if(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))>DSt){
							int MergedSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
							double *MergedPower=calloc(MergedSize,sizeof(double));
							int *MergedSE=calloc(MergedSize,sizeof(int));
							double *MergedDS=calloc(MergedSize,sizeof(double));
							int *MergedIndex;
							int ExceedUser=TotalServedMSs-ceil(Ntp*(1-BP));
							double ExceedPower=0.0;
							int ExceedDS=0;
							double ExceedSE=0.0;
							int m=0;
							if(ExceedUser>0){
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
										*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
										*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
										*(MergedDS+m)=*(DS+j+1);
										m++;
										}
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
										*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
										*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
										*(MergedDS+m)=*(DS+j+1);
										m++;
										}
								MergedIndex=Sq_sort(MergedDS,0,MergedSize-1,MergedSize);
								j=MergedSize-1;
								if(MergedSize>ExceedUser)
									while(j>MergedSize-ExceedUser-1){
										ExceedPower+=MergedPower[*(MergedIndex+j)];
										ExceedSE+=MergedSE[*(MergedIndex+j)];
										ExceedDS+=MergedDS[j];
										NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
										if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
											NewBS=*(SortedCumulativeObjectiveIndex+k);
											MaxIncrease=NewSE-obj;
											UserCutOff=2;
											j=-1;
										}
										else
											j--;
									
									}
								else if(MergedSize<=ExceedUser)
									while(j>-1){
										ExceedPower+=MergedPower[*(MergedIndex+j)];
										ExceedSE+=MergedSE[*(MergedIndex+j)];
										ExceedDS+=MergedDS[j];
										NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
										if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
											NewBS=*(SortedCumulativeObjectiveIndex+k);
											MaxIncrease=NewSE-obj;
											UserCutOff=2;
											j=-1;
										}
										else
											j--;
									
									}
							} 
							free(MergedPower);
							free(MergedSE);
							free(MergedDS);
// 							delete MergedIndex;
						}	
					}
					else if(NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))>MP+ABP){
						if(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))<=DSt){
							int MergedSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
							double *MergedPower=calloc(MergedSize,sizeof(double));
							int *MergedSE=calloc(MergedSize,sizeof(int));
							int *MergedIndex;
							int ExceedUser=TotalServedMSs-ceil(Ntp*(1-BP));
							double ExceedPower=0.0;
							double ExceedSE=0.0;
							int m=0;
							if(ExceedUser>0){
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
										*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
										*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
										m++;
										}
								for(j=0;j<Ntp;j++)
									if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
										*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
										*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
										m++;
										}
								MergedIndex=Sq_sort(MergedPower,0,MergedSize-1,MergedSize);
								j=MergedSize-1;
								if(MergedSize>ExceedUser)
									while(j>MergedSize-ExceedUser-1){
										ExceedPower+=MergedPower[j];
										ExceedSE+=MergedSE[*(MergedIndex+j)];
										NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
										if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP){
											NewBS=*(SortedCumulativeObjectiveIndex+k);
											MaxIncrease=NewSE-obj;
											UserCutOff=1;
											j=-1;
											printf("To New BS(%d) Transmission Power:%lf, Original BS(%d)\n",*(SortedCumulativeObjectiveIndex+NewBS)+1,NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower,*(SortedCumulativeObjectiveIndex+i)+1);
										}
										else
											j--;
									
									}
								else if(MergedSize<=ExceedUser)
									while(j>-1){
										ExceedPower+=MergedPower[j];
										ExceedSE+=MergedSE[*(MergedIndex+j)];
										NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
										if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP){
											NewBS=*(SortedCumulativeObjectiveIndex+k);
											MaxIncrease=NewSE-obj;
											UserCutOff=1;
											j=-1;
											printf("To New BS(%d) Transmission Power:%lf, Original BS(%d)\n",*(SortedCumulativeObjectiveIndex+NewBS)+1,NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower,*(SortedCumulativeObjectiveIndex+i)+1);
										}
										else
											j--;
									
									}
								if(j==MergedSize-ExceedUser){
									printf("Transition to New BS(%d) fail with Transmission Power:%lf and obj:%lf, Original BS(%d) obj:%lf\n",*(SortedCumulativeObjectiveIndex+k)+1,NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower,*(SortedCumulativeObjectiveIndex+i)+1,NewSE,obj);
								}
							} 
							free(MergedPower);
							free(MergedSE);
// 							delete MergedIndex;
						}
						else if(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))>DSt){
							double PR=(NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-MP)/MP;
							double DSR=(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-DSt)/DSt;
							if(PR>=DSR){
								int MergedSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
								double *MergedPower=calloc(MergedSize,sizeof(double));
								int *MergedSE=calloc(MergedSize,sizeof(int));
								double *MergedDS=calloc(MergedSize,sizeof(double));
								int *MergedIndex;
								int ExceedUser=TotalServedMSs-ceil(Ntp*(1-BP));
								double ExceedPower=0.0;
								int ExceedDS=0;
								double ExceedSE=0.0;
								int m=0;
								if(ExceedUser>0){
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
											*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
											*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
											*(MergedDS+m)=*(DS+j+1);
											m++;
											}
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
											*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
											*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
											*(MergedDS+m)=*(DS+j+1);
											m++;
											}
									MergedIndex=Sq_sort(MergedPower,0,MergedSize-1,MergedSize);
									j=MergedSize-1;
									if(MergedSize>ExceedUser)
										while(j>MergedSize-ExceedUser-1){
											ExceedPower+=MergedPower[j];
											ExceedSE+=MergedSE[*(MergedIndex+j)];
											ExceedDS+=MergedDS[*(MergedIndex+j)];
											NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
											if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
												NewBS=*(SortedCumulativeObjectiveIndex+k);
												MaxIncrease=NewSE-obj;
												UserCutOff=1;
												j=-1;
											}
											else
												j--;
										
										}
									else if(MergedSize<=ExceedUser)
										while(j>-1){
											ExceedPower+=MergedPower[j];
											ExceedSE+=MergedSE[*(MergedIndex+j)];
											ExceedDS+=MergedDS[*(MergedIndex+j)];
											NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
											if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
												NewBS=*(SortedCumulativeObjectiveIndex+k);
												MaxIncrease=NewSE-obj;
												UserCutOff=1;
												j=-1;
											}
											else
												j--;
										
										}
								} 
								free(MergedPower);
								free(MergedSE);
								free(MergedDS);
// 								delete MergedIndex;
							}
							else if(PR<DSR){
								int MergedSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
								double *MergedPower=calloc(MergedSize,sizeof(double));
								int *MergedSE=calloc(MergedSize,sizeof(int));
								double *MergedDS=calloc(MergedSize,sizeof(double));
								int *MergedIndex;
								int ExceedUser=TotalServedMSs-ceil(Ntp*(1-BP));
								double ExceedPower=0.0;
								int ExceedDS=0;
								double ExceedSE=0.0;
								int m=0;
								if(ExceedUser>0){
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
											*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
											*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
											*(MergedDS+m)=*(DS+j+1);
											m++;
											}
									for(j=0;j<Ntp;j++)
										if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
											*(MergedPower+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
											*(MergedSE+m)=(*(BR+j+1))/(*(DS+j+1));
											*(MergedDS+m)=*(DS+j+1);
											m++;
											}
									MergedIndex=Sq_sort(MergedDS,0,MergedSize-1,MergedSize);
									j=MergedSize-1;
									if(MergedSize>ExceedUser)
										while(j>MergedSize-ExceedUser-1){
											ExceedPower+=MergedPower[*(MergedIndex+j)];
											ExceedSE+=MergedSE[*(MergedIndex+j)];
											ExceedDS+=MergedDS[j];
											NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
											if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
												NewBS=*(SortedCumulativeObjectiveIndex+k);
												MaxIncrease=NewSE-obj;
												UserCutOff=2;
												j=-1;
											}
											else
												j--;
										
										}
									else if(MergedSize<=ExceedUser)
										while(j>-1){
											ExceedPower+=MergedPower[*(MergedIndex+j)];
											ExceedSE+=MergedSE[*(MergedIndex+j)];
											ExceedDS+=MergedDS[j];
											NewSE=(se-ExceedSE)/(p-*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+SBP+NewTransmissionPower-ExceedPower);
											if(NewSE-obj>MaxIncrease&&NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-ABP-ExceedPower<=MP&&*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))+*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k))-ExceedDS<=DSt){
												NewBS=*(SortedCumulativeObjectiveIndex+k);
												MaxIncrease=NewSE-obj;
												UserCutOff=2;
												j=-1;
											}
											else
												j--;
										
										}
								} 
								free(MergedPower);
								free(MergedSE);
								free(MergedDS);
// 								delete MergedIndex;
							}
						}
					}
				}/*
				if(((*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i)))-NewObjective)<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+SBP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					printf("1. Swith to New BS:%d\n",NewBS+1);
					printf("1. MaxObjective:%lf\n",MaxNewObjective);
					UserCutOff=0;
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))>0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					printf("2. Swith to New BS:%d\n",NewBS+1);
					printf("2. MaxObjective:%lf\n",MaxNewObjective);
					UserCutOff=0;
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					printf("3. Swith to New BS:%d\n",NewBS+1);
					printf("3. MaxObjective:%lf\n",MaxNewObjective);
					UserCutOff=0;
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))>MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					printf("4. Try to switch to New BS:%d....\n",*(SortedCumulativeObjectiveIndex+k)+1);
					int IntegralBSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
					double IntegralBS[IntegralBSize];//=calloc(IntegralBSize,sizeof(double));
			 		double IntegralObjective[IntegralBSize];//=calloc(IntegralBSize,sizeof(double));
					ArrayInitialization(IntegralBS,IntegralBSize);
			 		ArrayInitialization(IntegralObjective,IntegralBSize);
			 		int *SortedIntegralBSPowerIndex;
					int BlockableUserNumber=TotalServedMSs-ceil(Ntp*(1-BP));
					double SupressPower=0.0;
			 		double CalculatedObjective=0.0;
					int m=0;
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
// 							for(m=0;m<*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m++){
								*(IntegralBS+m)=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
								printf("%d: %lf  ",m,*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1));
							m++;
							}
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
							//for(m=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m<IntegralBSize;m++){
								*(IntegralBS+m)=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
							m++;
							}
					
					SortedIntegralBSPowerIndex=Sq_sort(IntegralBS,0,IntegralBSize-1,IntegralBSize);
					printf("%d\n",IntegralBSize);
					
					for(j=IntegralBSize-1;j>IntegralBSize-BlockableUserNumber;j--){
						SupressPower+=*(IntegralBS+j);
						CalculatedObjective+=IntegralObjective[*(SortedIntegralBSPowerIndex+j)];
					}
					
					printf("suppress %lf, %lf\n",SupressPower,NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-SupressPower);
					//free(IntegralObjective);
					if(NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-SupressPower<=MP+ABP){
						NewBS=*(SortedCumulativeObjectiveIndex+k);
						MaxNewObjective=NewObjective-CalculatedObjective;
						UserCutOff=1;
						NewPower=NewTransmissionPower;
						printf("4. DO Swith to New BS:%d, NewBS:%d\n",*(SortedCumulativeObjectiveIndex+k)+1,NewBS);
					}
					printf("4. MaxObjective:%lf\n",MaxNewObjective);
					//free(IntegralBS);
					
					printf("NewBS:%d, MaxNewObjective:%lf, UserCutOff:%d, NewPower:%lf\n",NewBS,MaxNewObjective,UserCutOff,NewPower);
					//printf("suppress %lf, %lf\n",SupressPower,NewPower+*(UnsortedCumulativePower+NewBS)-SupressPower);
					
				}
				else if(((*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i)))-NewObjective)<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))>MP+SBP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0)){
					printf("5. Consider switch to New BS:%d....\n",*(SortedCumulativeObjectiveIndex+k)+1);
					int IntegralBSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
					double *IntegralBS=malloc(IntegralBSize*sizeof(double));
			 		double *IntegralObjective=malloc(IntegralBSize*sizeof(double));
					ArrayInitialization(IntegralBS,IntegralBSize);
			 		ArrayInitialization(IntegralObjective,IntegralBSize);
			 		int *SortedIntegralBSPowerIndex;
					int BlockableUserNumber=TotalServedMSs-ceil(Ntp*(1-BP));
					double SupressPower=0.0;
			 		double CalculatedObjective=0.0;
					int m=0;;
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
							//for(m=0;m<*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m++){
								IntegralBS[m]=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
							m++;
							}
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
// 							for(m=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m<IntegralBSize;m++){
								IntegralBS[m]=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
							m++;
							}
					SortedIntegralBSPowerIndex=Sq_sort(IntegralBS,0,IntegralBSize-1,IntegralBSize);
					for(j=IntegralBSize-1;j>IntegralBSize-BlockableUserNumber;j--){
						SupressPower+=IntegralBS[j];
						CalculatedObjective+=IntegralObjective[*(SortedIntegralBSPowerIndex+j)];
					}
					free(IntegralBS);
					free(IntegralObjective);
					if(NewTransmissionPower-SupressPower<=MP&&NewObjective-CalculatedObjective>MaxNewObjective){
						NewBS=*(SortedCumulativeObjectiveIndex+k);
						MaxNewObjective=NewObjective-CalculatedObjective;
						UserCutOff=1;
						NewPower=NewTransmissionPower;
						printf("5. DO Swith to New BS:%d, NewBS:%d\n",*(SortedCumulativeObjectiveIndex+k)+1,NewBS);
					}
					printf("5. MaxObjective:%lf\n",MaxNewObjective);
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))>0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))>MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					printf("6. Swith to New BS:%d\n",*(SortedCumulativeObjectiveIndex+k)+1);
					int IntegralBSize=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+k))+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
					double *IntegralBS=malloc(IntegralBSize*sizeof(double));
			 		double *IntegralObjective=malloc(IntegralBSize*sizeof(double));
					ArrayInitialization(IntegralBS,IntegralBSize);
			 		ArrayInitialization(IntegralObjective,IntegralBSize);
			 		int *SortedIntegralBSPowerIndex;
					int BlockableUserNumber=TotalServedMSs-ceil(Ntp*(1-BP));
					double SupressPower=0.0;
			 		double CalculatedObjective=0.0;
					int m=0;
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
// 							for(m=0;m<*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m++){
								IntegralBS[m]=*(P+((*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1))+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
							m++;
							}
					for(j=0;j<Ntp;j++)
						if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+k)){
// 							for(m=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));m<IntegralBSize;m++){
								IntegralBS[m]=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
			 					IntegralObjective[m]=UnsortedObjective[*(SortedCumulativeObjectiveIndex+k)][j];
							m++;
							}
					SortedIntegralBSPowerIndex=Sq_sort(IntegralBS,0,IntegralBSize-1,IntegralBSize);
					for(j=IntegralBSize-1;j>IntegralBSize-BlockableUserNumber;j--){
						SupressPower+=IntegralBS[j];
						CalculatedObjective+=IntegralObjective[*(SortedIntegralBSPowerIndex+j)];
					}
					free(IntegralBS);
					free(IntegralObjective);
					if(NewTransmissionPower+*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))-SupressPower<=MP+ABP){
						NewBS=*(SortedCumulativeObjectiveIndex+k);
						MaxNewObjective=NewObjective-CalculatedObjective;
						UserCutOff=1;
						NewPower=NewTransmissionPower;
						printf("6. DO Swith to New BS:%d, NewBS:%d\n",*(SortedCumulativeObjectiveIndex+k)+1,NewBS);
					}
					printf("6. MaxObjective:%lf\n",MaxNewObjective);
				}*/
			}		
			//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
			//printf("MaxObjective:%lf\n",MaxNewObjective);
		}	
	}
// 	for(l=0;l<Tbs;l++)
// 		printf("BS%d is in status %d.\n",l+1,BStatus[l]);
	//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
	if(NewBS!=-1&&UserCutOff==0){
		double Energy;
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i))
				Energy+=*(P+((NewBS+1)*(Ntp+1))+j+1);
		printf("Users in BS%d(%d) is gonna be reassigned to BS%d(%d) with an additional power level of %lf\n",*(SortedCumulativeObjectiveIndex+i)+1,*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i)),NewBS+1,*(CumulativeServingMSs+NewBS),Energy);
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){	
				if(*(CumulativeServingMSs+NewBS)==0)
					*(UnsortedCumulativePower+NewBS)+=ABP-SBP;						
				*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
				*(UnsortedCumulativeObjective+NewBS)+=UnsortedObjective[NewBS][j];
				*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
				*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=UnsortedObjective[*(SortedCumulativeObjectiveIndex+i)][j];
				*(CumulativeDS+NewBS)+=*(DS+j+1);
				*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
				*(CumulativeServingMSs+NewBS)+=1;
				*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
				*(temp_dBM+j)=NewBS;
				*(CumulativeBR+NewBS)+=*(BR+j+1);
				*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);						
				//if(BStatus)
										

			}
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		if(*(CumulativeServingMSs+NewBS)==0||*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))!=0){
			printf("j. Program Stops! BS%d serves for wrong number of users!\n",*(SortedCumulativeObjectiveIndex+i));
			exit(1);
		}
		BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
		*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))=SBP;
		if((*(UnsortedCumulativePower+NewBS)-ABP)==MP||(*(CumulativeDS+NewBS)==DSt))
			BStatus[NewBS]=2;
		else
			BStatus[NewBS]=1;
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		/*for(l=0;l<Tbs;l++)
			printf("BS%d is in status %d.\n",l+1,BStatus[l]);*/
		i=0;
		for(l=0;l<Tbs;l++)
			if(BStatus[l]==0){
				if(*(UnsortedCumulativePower+l)-SBP>10E-8){
					printf("i1. rogram Stops! BS%d(%d) transmits power at %g level with DS%d %lf.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l),*(UnsortedCumulativePower+l)-SBP);
					exit(1);
				}
				else if(fabs(*(CumulativeDS+l)-0)>0.0000000000000000000000000001){
					printf("i2. rogram Stops! BS%d(%d) transmits power at %g level with DS%d %lf.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l),*(UnsortedCumulativePower+l)-SBP);
					exit(1);
				}
			}
			else 
			if(BStatus[l]!=0&&(*(UnsortedCumulativePower+l)>MP+ABP||*(CumulativeDS+l)>DSt)){
				printf("h. Program Stops! BS%d(%d) transmits power at %lf level with DS%d.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
		SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
	}
	else 
	if(NewBS!=-1&&UserCutOff==1){
// 		printf("The program gets here.\n");
		/*int IntegralBSize=*(CumulativeServingMSs+NewBS)+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
		double Integral[IntegralBSize];//=calloc(IntegralBSize,sizeof(double));
// 		double *IntegralObjective=malloc(IntegralBSize*sizeof(double));
		ArrayInitialization(Integral,IntegralBSize);
// 		ArrayInitialization(IntegralObjective,IntegralBSize);
// 		int *SortedIntegralBSPowerIndex;
		int BlockableUserNumber=TotalServedMSs-ceil(Ntp*(1-BP));
		double SupressPower=0.0;
// 		double CalculatedObjective=0.0;
 		int k=0;
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
// 				for(k=0;k<*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));k++){
					Integral[k]=*(P+((NewBS+1)*(Ntp+1))+j+1);
 					printf("%d. %lf\n",k,Integral[k]);
// 					IntegralObjective[k]=UnsortedObjective[NewBS][j];
				k++;
				}
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==NewBS){
// 				for(k=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));k<IntegralBSize;k++){
					Integral[k]=*(P+(NewBS+1)*(Ntp+1)+j+1);
// 					IntegralObjective[k]=UnsortedObjective[NewBS][j];
				k++;
				}
		//SortedIntegralBSPowerIndex=
		Sq_sort(Integral,0,IntegralBSize-1,IntegralBSize);
		for(j=IntegralBSize-1;j>IntegralBSize-BlockableUserNumber;j--){
			SupressPower+=Integral[j];
			//CalculatedObjective+=InegralObjective[*(SortedIntegralBSPowerIndex+j)];
		}printf("1. The program gets here.\n");
		//free(Integral);
		//free(IntegralObjective);
		printf("suppress %lf, %lf\n",SupressPower,NewPower+*(UnsortedCumulativePower+NewBS)-SupressPower);*/
// 		if((NewPower+*(UnsortedCumulativePower+NewBS)-SupressPower<=MP+ABP&&BStatus[NewBS]==1)||(NewPower-SupressPower<=MP&&BStatus[NewBS]==0)){
		if(MaxIncrease>0){
			for(j=0;j<Ntp;j++)
				if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
					if(*(CumulativeServingMSs+NewBS)==0)
						*(UnsortedCumulativePower+NewBS)+=ABP-SBP;
					*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
					*(UnsortedCumulativeObjective+NewBS)+=UnsortedObjective[NewBS][j];
					
					
					*(CumulativeDS+NewBS)+=*(DS+j+1);
					
					*(CumulativeServingMSs+NewBS)+=1;
					
					*(temp_dBM+j)=NewBS;
					*(CumulativeBR+NewBS)+=*(BR+j+1);
					
					
						*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
						*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=UnsortedObjective[*(SortedCumulativeObjectiveIndex+i)][j];
						*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
						*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
						*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);	
					
				}
			if(*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))==0){
				BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
				*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))=SBP;
			}
			else{
				printf("g. Program Stops! BS%d transmits power at %lf level.\n",*(SortedCumulativeObjectiveIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i)));
				exit(1);
			}
			IndexMaxP=-1;
			MaxP=0.0;
			while(*(UnsortedCumulativePower+NewBS)>MP+ABP||*(CumulativeDS+NewBS)>DSt){
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==NewBS&&*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)>MaxP){
						MaxP=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1);
						IndexMaxP=j;
					}
					else
					if(*(temp_dBM+j)==NewBS&&*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)==MaxP)
						if(UnsortedObjective[NewBS][IndexMaxP]>UnsortedObjective[NewBS][j]){
							if(MaxP!=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)){
								printf("3. Program Stops!\n");
								exit(1);
							}
							MaxP=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1);
							IndexMaxP=j;
					}
				}
				if(IndexMaxP!=-1){
					*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
					*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
					CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
					TotalServedMSs--;
					*(CumulativeBR+*(temp_dBM+IndexMaxP))-=*(BR+IndexMaxP+1);
					*(UnsortedCumulativeObjective+*(temp_dBM+IndexMaxP))-=UnsortedObjective[*(temp_dBM+IndexMaxP)][IndexMaxP];
					if(*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))==MP+ABP||(*(CumulativeDS+*(temp_dBM+IndexMaxP))==DSt))
						BStatus[*(temp_dBM+IndexMaxP)]=2;
					else
						BStatus[*(temp_dBM+IndexMaxP)]=1;
					*(temp_dBM+IndexMaxP)=-1;
				}
				IndexMaxP=-1;
				MaxP=0.0;
			}	
			
		i=0;
		for(l=0;l<Tbs;l++){
			if(BStatus[l]==0){
				if(*(UnsortedCumulativePower+l)-SBP>10E-8){
					printf("f. Program Stops! BS%d(%d) transmits power at %lf level with DS%d.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l));
					exit(1);
				}
				else if(*(CumulativeDS+l)!=0){
				printf("e. Program Stops! BS%d(%d) transmits power at %lf level with DS%d.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l));
				exit(1);
				}
			}
			else 
			if(BStatus[l]!=0&&(*(UnsortedCumulativePower+l)>MP+ABP||*(CumulativeDS+l)>DSt)){
				printf("d. Program Stops! BS%d transmits power at %lf level with DS%d.\n",l+1,*(UnsortedCumulativePower+l),*(CumulativeDS+l));
				exit(1);
			}
		}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
		SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
		}
		UserCutOff=0;
		
	}
	else 
	if(NewBS!=-1&&UserCutOff==2){
		if(MaxIncrease>0){
			for(j=0;j<Ntp;j++)
				if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
					*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
					*(UnsortedCumulativeObjective+NewBS)+=UnsortedObjective[NewBS][j];
					*(CumulativeDS+NewBS)+=*(DS+j+1);					
					*(CumulativeServingMSs+NewBS)+=1;
					*(temp_dBM+j)=NewBS;
					*(CumulativeBR+NewBS)+=*(BR+j+1);
					*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
					*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=UnsortedObjective[*(SortedCumulativeObjectiveIndex+i)][j];
					*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
					*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
					*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);		
				}
			if(*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))==0){
				BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
				*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))=SBP;
			}
			else{
				printf("Program Stops! BS%d transmits power at %lf level.\n",*(SortedCumulativeObjectiveIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i)));
				exit(1);
			}
			int IndexMaxDSs=-1;
			int MaxDSs=0;
			while(*(UnsortedCumulativePower+NewBS)>MP+ABP||*(CumulativeDS+NewBS)>DSt){
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==NewBS&&*(DS+j+1)>MaxDSs){
						MaxDSs=*(DS+j+1);
						IndexMaxDSs=j;
					}
					else
					if(*(temp_dBM+j)==NewBS&&*(DS+j+1)==MaxDSs)
						if(UnsortedObjective[NewBS][IndexMaxP]>UnsortedObjective[NewBS][j]){
							if(MaxDSs!=*(DS+j+1)){
								printf("3. Program Stops!\n");
								exit(1);
							}
							MaxDSs=*(DS+j+1);
							IndexMaxDSs=j;
					}
				}
				if(IndexMaxDSs!=-1){
					*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
					*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
					CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
					TotalServedMSs--;
					*(CumulativeBR+*(temp_dBM+IndexMaxP))-=*(BR+IndexMaxP+1);
					*(UnsortedCumulativeObjective+*(temp_dBM+IndexMaxP))-=UnsortedObjective[*(temp_dBM+IndexMaxP)][IndexMaxP];
					if(*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))==MP+ABP||(*(CumulativeDS+*(temp_dBM+IndexMaxP))==DSt))
						BStatus[*(temp_dBM+IndexMaxP)]=2;
					else
						BStatus[*(temp_dBM+IndexMaxP)]=1;
					*(temp_dBM+IndexMaxP)=-1;
				}
				IndexMaxP=-1;
				MaxP=0.0;
			}	
			
		i=0;
		for(l=0;l<Tbs;l++){
			if(BStatus[l]==0){
				if(*(UnsortedCumulativePower+l)!=SBP){
					printf("c. Program Stops! BS%d(%d) transmits power at %lf level with DS%d.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l));
					exit(1);
				}
				else if(*(CumulativeDS+l)!=0){
				printf("b. Program Stops! BS%d(%d) transmits power at %lf level with DS%d.\n",l+1,BStatus[l],*(UnsortedCumulativePower+l),*(CumulativeDS+l));
				exit(1);
				}
			}
			else 
			if(BStatus[l]!=0&&(*(UnsortedCumulativePower+l)>MP+ABP||*(CumulativeDS+l)>DSt)){
				printf("a. Program Stops! BS%d transmits power at %lf level with DS%d.\n",l+1,*(UnsortedCumulativePower+l),*(CumulativeDS+l));
				exit(1);
			}
		}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
		SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
		}
	}		
		
	/*for(l=0;l<Tbs;l++)
		printf("%d ",*(SortedCumulativePowerIndex+l)+1);
		printf("\n");
	for(l=0;l<Tbs;l++)
		printf("BS%d is in status %d.\n",l+1,BStatus[l]);	*/	
	
}/*
double NewObjective=0.0;
double NewTransmissionPower=0.0;
double NewPower=0.0;
double MaxNewObjective;
int NewBS=-1;
int UserCutOff=0;
double MaxP=0.0;
int IndexMaxP=-1;
for(i=0;i<Tbs-1;i++){
	NewBS=-1;
	MaxNewObjective=0.0;
	NewPower=0.0;
	if(BStatus[*(SortedCumulativeObjectiveIndex+i)]==1||BStatus[*(SortedCumulativeObjectiveIndex+i)]==2){
		//printf("Old BS:%d\n",*(SortedCumulativePowerIndex+i));
		for(k=Tbs-1;k>-1;k--){
			if(*(SortedCumulativeObjectiveIndex+k)!=*(SortedCumulativeObjectiveIndex+i)){
				NewTransmissionPower=0.0;
				NewObjective=0.0;
				//printf("Old BS:%d, New BS:%d\n",*(SortedCumulativePowerIndex+i)+1,*(SortedCumulativePowerIndex+k)+1);
				if((BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)||(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0))
					if(((*(CumulativeDS+*(SortedCumulativeObjectiveIndex+k)))+(*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))))>DSt)
						for(j=0;j<Ntp;j++)
							if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
								NewTransmissionPower+=*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1);
								NewObjective+=(*(BR+j+1))/(*(DS+j+1))/(*(P+(*(SortedCumulativeObjectiveIndex+k)+1)*(Ntp+1)+j+1));
							}
				//printf("To New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",*(SortedCumulativePowerIndex+k)+1,NewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));		
				if(((*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i)))-NewObjective)<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+SBP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					//printf("Swith to New BS:%d\n",NewBS);
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))>0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					//printf("Swith to New BS:%d\n",NewBS);
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					//printf("Swith to New BS:%d\n",NewBS);
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))>MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					UserCutOff=1;
					NewPower=NewTransmissionPower;
					//printf("Swith to New BS:%d\n",NewBS);
				}
				else if(((*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i)))-NewObjective)<0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))>MP+SBP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==0)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					UserCutOff=1;
					NewPower=NewTransmissionPower;
					//printf("Swith to New BS:%d\n",NewBS);
				}
				else if((NewObjective-(*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))))>0&&MaxNewObjective<NewObjective&&(NewTransmissionPower+(*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+k))))<=MP+ABP&&(BStatus[*(SortedCumulativeObjectiveIndex+k)]==1)){
					NewBS=*(SortedCumulativeObjectiveIndex+k);
					MaxNewObjective=NewObjective;
					UserCutOff=1;
					NewPower=NewTransmissionPower;
			}		
			//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
			
		}	
	}
// 	for(l=0;l<Tbs;l++)
// 		printf("BS%d is in status %d.\n",l+1,BStatus[l]);
	//printf("New BS(%d) Transmission Power:%lf, Original BS(%d) Transmission Power:%lf\n",NewBS+1,MinNewTransmissionPower,*(SortedCumulativePowerIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativePowerIndex+i)));
	if(NewBS!=-1&&UserCutOff==0){
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){	
				if(*(CumulativeServingMSs+NewBS)==0)
					*(UnsortedCumulativePower+NewBS)+=ABP-SBP;						
				*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
				*(UnsortedCumulativeObjective+NewBS)+=UnsortedObjective[NewBS][j];
				*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
				*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=UnsortedObjective[*(SortedCumulativeObjectiveIndex+i)][j];
				*(CumulativeDS+NewBS)+=*(DS+j+1);
				*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
				*(CumulativeServingMSs+NewBS)+=1;
				*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
				*(temp_dBM+j)=NewBS;
				*(CumulativeBR+NewBS)+=*(BR+j+1);
				*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);						

										

			}
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		if(*(CumulativeServingMSs+NewBS)==0||*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))!=0){
			printf("Program Stops! BS%d serves for wrong number of users!\n",*(SortedCumulativeObjectiveIndex+i));
			exit(1);
		}
		BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
		*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+=SBP-ABP;
		if((*(UnsortedCumulativePower+NewBS)-ABP)==MP||(*(CumulativeDS+NewBS)==DSt))
			BStatus[NewBS]=2;
		else
			BStatus[NewBS]=1;
		//printf("BS%d(%d) Swith to New BS:%d(%d)\n",*(SortedCumulativePowerIndex+i)+1,BStatus[*(SortedCumulativePowerIndex+i)],NewBS+1,BStatus[NewBS]);
		
		i=0;
		for(l=0;l<Tbs;l++)
			if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
				printf("Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
		SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
	}
	else 
	if(NewBS!=-1&&UserCutOff==1){
		int IntegralBSize=*(CumulativeServingMSs+NewBS)+*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));
		double *IntegralBS=malloc(IntegralBSize*sizeof(double));
		
		ArrayInitialization(IntegralBS,IntegralBSize);
		//int *SortedIntegralBSPowerIndex;
		int BlockableUserNumber=TotalServedMSs-ceil(Ntp*(1-BP));
		double SupressPower=0.0;
		
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i))
				for(k=0;k<*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));k++)
					IntegralBS[k]=*(P+((NewBS+1)*(Ntp+1))+j+1);
		for(j=0;j<Ntp;j++)
			if(*(temp_dBM+j)==NewBS)
				for(k=*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i));k<IntegralBSize;k++)
					IntegralBS[k]=*(P+(NewBS+1)*(Ntp+1)+j+1);
		//SortedIntegralBSPowerIndex=
		Sq_sort(IntegralBS,0,IntegralBSize-1,IntegralBSize);
		for(j=IntegralBSize-1;j>IntegralBSize-BlockableUserNumber;j--)
			SupressPower+=IntegralBS[j];
		free(IntegralBS);
		if(NewPower+*(UnsortedCumulativePower+NewBS)-SupressPower<=MP){
			for(j=0;j<Ntp;j++)
				if(*(temp_dBM+j)==*(SortedCumulativeObjectiveIndex+i)){
					*(UnsortedCumulativePower+NewBS)+=*(P+((NewBS+1)*(Ntp+1))+j+1);
					*(UnsortedCumulativeObjective+NewBS)+=UnsortedObjective[NewBS][j];
					*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))-=*(P+(*(SortedCumulativeObjectiveIndex+i)+1)*(Ntp+1)+j+1);
					*(UnsortedCumulativeObjective+*(SortedCumulativeObjectiveIndex+i))-=UnsortedObjective[*(SortedCumulativeObjectiveIndex+i)][j];
					*(CumulativeDS+NewBS)+=*(DS+j+1);
					*(CumulativeDS+*(SortedCumulativeObjectiveIndex+i))-=*(DS+j+1);
					*(CumulativeServingMSs+NewBS)+=1;
					*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))-=1;
					*(temp_dBM+j)=NewBS;
					*(CumulativeBR+NewBS)+=*(BR+j+1);
					*(CumulativeBR+*(SortedCumulativeObjectiveIndex+i))-=*(BR+j+1);	
				}
			if(*(CumulativeServingMSs+*(SortedCumulativeObjectiveIndex+i))==0){
				BStatus[*(SortedCumulativeObjectiveIndex+i)]=0;
				*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i))+=SBP-ABP;
			}
			else{
				printf("Program Stops! BS%d transmits power at %lf level.\n",*(SortedCumulativeObjectiveIndex+i)+1,*(UnsortedCumulativePower+*(SortedCumulativeObjectiveIndex+i)));
				exit(1);
			}
			IndexMaxP=-1;
			MaxP=0.0;
			while(*(UnsortedCumulativePower+NewBS)>MP){
				for(j=0;j<Ntp;j++){
					if(*(temp_dBM+j)==NewBS&&*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)>MaxP){
						MaxP=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1);
						IndexMaxP=j;
					}
					else
					if(*(temp_dBM+j)==NewBS&&*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)==MaxP)
						if(UnsortedObjective[NewBS][IndexMaxP]>UnsortedObjective[NewBS][j])){
							if(MaxP!=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1)){
								printf("Program Stops!\n");
								exit(1);
							}
							MaxP=*(P+(*(temp_dBM+j)+1)*(Ntp+1)+j+1);
							IndexMaxP=j;
					}
				}
				if(IndexMaxP!=-1){
					*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))-=*(P+(*(temp_dBM+IndexMaxP)+1)*(Ntp+1)+IndexMaxP+1);
					*(CumulativeDS+*(temp_dBM+IndexMaxP))-=*(DS+IndexMaxP+1);
					CumulativeServingMSs[*(temp_dBM+IndexMaxP)]--;
					TotalServedMSs--;
					*(CumulativeBR+*(temp_dBM+IndexMaxP))-=*(BR+IndexMaxP+1);
					*(UnsortedCumulativeObjective+*(temp_dBM+IndexMaxP))-=UnsortedObjective[*(temp_dBM+IndexMaxP)][IndexMaxP];
					if(*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP))==MP+ABP||(*(CumulativeDS+*(temp_dBM+IndexMaxP))==DSt))
						BStatus[*(temp_dBM+IndexMaxP)]=2;
					else
						BStatus[*(temp_dBM+IndexMaxP)]=1;
					*(temp_dBM+IndexMaxP)=-1;
				}
				IndexMaxP=-1;
				MaxP=0;
			}	
			
		i=0;
		for(l=0;l<Tbs;l++)
			if(*(UnsortedCumulativePower+l)>MP||*(CumulativeDS+l)>DSt){
				printf("Program Stops! BS%d transmits power at %lf level.\n",l+1,*(UnsortedCumulativePower+l));
				exit(1);
			}
		for(l=0;l<Tbs;l++)
			*(SortedCumulativeObjective+l)=*(UnsortedCumulativeObjective+l);
		SortedCumulativeObjectiveIndex=Sq_sort(SortedCumulativeObjective,0,Tbs-1,Tbs);
		}
		UserCutOff=0;
		
	}
			
		
		
	
}	*//*			
while(TotalServedMSs>ceil(Ntp*(1-BP))){
	//for(i=0;i<Ntp;i++)		printf("MS%d, index: %d, power: %lf\n ",i+1,*(temp_dBM+i)+1,*(P+(*(temp_dBM+i)+1)*(Ntp+1)+i+1));
	for(i=0;i<Ntp;i++){
		if(*(temp_dBM+i)!=-1&&UnsortedObjective[*(temp_dBM+i)][i]<Min){
			Min=UnsortedObjective[*(temp_dBM+i)][i];
			IndexMin=i;
		}
	}
	if(IndexMin!=-1){
		*(UnsortedCumulativePower+*(temp_dBM+IndexMin))-=*(P+(*(temp_dBM+IndexMin)+1)*(Ntp+1)+IndexMin+1);
		*(UnsortedCumulativeObjective+*(temp_dBM+IndexMin))-=UnsortedObjective[*(temp_dBM+IndexMin)][IndexMin];
		*(CumulativeDS+*(temp_dBM+IndexMin))-=*(DS+IndexMin+1);
		*(CumulativeBR+*(temp_dBM+IndexMin))-=*(BR+IndexMin+1);
		CumulativeServingMSs[*(temp_dBM+IndexMin)]--;
		TotalServedMSs--;
//printf("BS %d SERVES # %d\n",*(temp_dBM+IndexMaxP),CumulativeServingMSs[*(temp_dBM+IndexMaxP)]);
		
		if(CumulativeServingMSs[*(temp_dBM+IndexMin)]==0)
			BStatus[*(temp_dBM+IndexMin)]=0;
		else
			BStatus[*(temp_dBM+IndexMin)]=1;
		*(temp_dBM+IndexMin)=-1;
	}
	IndexMin=-1;
	Min=1000.0;
	//printf("REMOVE INDEX%d, %lf",IndexMaxP,*(UnsortedCumulativePower+*(temp_dBM+IndexMaxP)));
}*/
temp=0;
for(i=0;i<Tbs;i++)
	temp+=CumulativeServingMSs[i];
if(TotalServedMSs!=temp){
printf("4. Inconsistent numbers!%d %d",temp,TotalServedMSs);
exit(1);}

for(j=0;j<Tbs;j++){
	if((*(UnsortedCumulativePower+j)>(MP+ABP))||(*(CumulativeDS+j)>DSt)){
		printf("Constraint violation! Program Stops!!\n");
		printf("BS%d Power %lf DS %d\n",j,*(UnsortedCumulativePower+j),*(CumulativeDS+j)>DSt);
		exit(1);
	}
}

	for(i=0;i<Ntp;i++)
		if(*(temp_dBM+i)!=-1)
			dBM[*(temp_dBM+i)][i]=1;
temp=0;
for(i=0;i<Tbs;i++)
	temp+=CumulativeServingMSs[i];
if(TotalServedMSs!=temp){
printf("6. Inconsistent numbers!%d %d",temp,TotalServedMSs);
exit(1);}

       

if(TotalServedMSs<ceil(Ntp*(1-BP))){
	printf("error!!Possible Served MSs:%d. ", TotalServedMSs);
	printf("The required number(%g) of total served MSs(%d) can not be achieved!\n", ceil(Ntp*(1-BP)),Ntp);
	exit(1);
}
       /*
       for(i=0;i<Tbs;i++){
       *(SortedCumulativePower+i)=*(UnsortedCumulativePower+i);

       }
       printf("\n");
       
       SortedCumulativePowerIndex=q_sort(SortedCumulativePower,0,Tbs-1,Tbs);
       for(i=0;i<Tbs;i++)
       printf("Unsorted Cumulative Power of BS %d: %lf \n",i+1,*(UnsortedCumulativePower+i));
       printf("\n");
       for(i=0;i<Tbs;i++)
       printf("%d. Sorted Cumulative Power of BS %d: %lf \n",i+1,*(SortedCumulativePowerIndex+i)+1,*(SortedCumulativePower+i));
       printf("\n");
       for(i=0;i<Tbs;i++)
       printf("%d ",*(SortedCumulativePowerIndex+i)+1);
       printf("\n");
       
       for(k=Tbs-1;k>-1;k--){
                           if(TotalServedMSs-*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))>ceil(Ntp*(1-BP))){
                           for(i=0;i<Ntp;i++)
                           *(*dBM+(*(SortedCumulativePowerIndex+k)*Ntp)+i)=0;
                           *(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))=0;
                           *(CumulativeDS+*(SortedCumulativePowerIndex+k))=0;
			   TotalServedMSs-=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k));
*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))=0;
                           
//printf("Total served MSs: %d\n",TotalServedMSs);
                           }
                           else if(TotalServedMSs-*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))==ceil(Ntp*(1-BP))){
                                for(i=0;i<Ntp;i++)
                                *(*dBM+(*(SortedCumulativePowerIndex+k)*Ntp)+i)=0;
                                *(UnsortedCumulativePower+*(SortedCumulativePowerIndex+k))=0;
                                *(CumulativeDS+*(SortedCumulativePowerIndex+k))=0;
				TotalServedMSs-=*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k));
*(CumulativeServingMSs+*(SortedCumulativePowerIndex+k))=0;
                                
//printf("Total served MSs: %d\n",TotalServedMSs);
                                break;
                                }
                           }*/
                           

for(i=0;i<Tbs;i++){
       printf("Power from BS %d: %lf\n",i+1,*(UnsortedCumulativePower+i));
        TotalPower+=*(UnsortedCumulativePower+i);
       printf("DS from BS %d: %d\n",i+1,*(CumulativeDS+i));
	printf("BR from BS %d: %lf\n",i+1,*(CumulativeBR+i));
	printf("Objective from BS %d: %lf\n",i+1,*(UnsortedCumulativeObjective+i));
      TotalDSs+=*(CumulativeDS+i);
 	TotalBR+=*(CumulativeBR+i);
	ObjectiveValue+=*(UnsortedCumulativeObjective+i);
       printf("MSs served by BS %d: %d\n",i+1,*(CumulativeServingMSs+i));
       for(j=0;j<Ntp;j++){
       printf("dBM%d_%d: %d\n",i+1,j+1,dBM[i][j]);
		if(dBM[i][j]==1)
			NOV+=(BR[j+1]/DS[j+1]);
	}
       }
       
       printf("Total Power %lf\n",TotalPower);
       printf("Total DSs: %d\n",TotalDSs);
	printf("Total BR: %lf\n",TotalBR);
       printf("Total served MSs: %d\n",TotalServedMSs);
	printf("Objective Value: %lf\n",ObjectiveValue);
       for(j=1;j<=Ntp;j++)
       printf("MS %d needs DSs: %d\n",j,*(DS+j));
       
       SdBMprinter(dBM[0],Ntp,Tbs,DS,P);
       for(i=0;i<Tbs;i++){
       printf("Power from BS %d: %lf\n",i+1,*(UnsortedCumulativePower+i));
       printf("DS from BS %d: %d\n",i+1,*(CumulativeDS+i));
	printf("BR from BS %d: %lf\n",i+1,*(CumulativeBR+i));
       printf("MSs served by BS %d: %d\n",i+1,*(CumulativeServingMSs+i)); 
	printf("Objective from BS %d: %lf\n\n",i+1,*(UnsortedCumulativeObjective+i));
       }
       printf("Total Power %lf\n",TotalPower);
       printf("Total DSs: %d\n",TotalDSs);
	printf("Total BR: %lf\n",TotalBR);
       printf("Total served MSs: %d\n",TotalServedMSs);
	printf("Objective Value: %lf\n",TotalBR-TotalDSs-TotalPower);
	printf("New Objective Value: %lf\n",NOV/TotalPower);
       
       end_tick=clock();
       elapsed=(double)(end_tick-start_tick)/CLOCKS_PER_SEC;
       printf("Running Time: %.9f",elapsed);

	coverage=SCoverageFinder(dBM[0],ptrD,Tbs,Ntp);

FILE *HeuristicOut;
if ((HeuristicOut=fopen("SHeuristicOut.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen SHeuristicOut.txt successfully!\n");

fprintf(HeuristicOut,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
for(i=0;i<Tbs;i++)
fprintf(HeuristicOut," \"BS[%d]\" %d %d %g\n",i+1,Xbs[i+1],Ybs[i+1],*(coverage+i));
fprintf(HeuristicOut,"\n\n");
               for(i=0;i<Tbs;i++){
fprintf(HeuristicOut,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
               for(j=0;j<Ntp;j++)
		if(dBM[i][j]==1){
fprintf(HeuristicOut,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
fprintf(HeuristicOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW\n",i+1,Xbs[i+1],Ybs[i+1],j+1,Xtp[j+1],Ytp[j+1],DS[j+1],*(P+(i+1)*(Ntp+1)+j+1));
}
fprintf(HeuristicOut,"\n\n");
}
fprintf(HeuristicOut," plot \"SHeuristicOut.txt\" index 0:0 using 2:3:1 notitle with labels, \"SHeuristicOut.txt\" index 0:0 using 2:3:4 title \"cell size\" with circles, \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
for(i=1;i<=Tbs;i++)
fprintf(HeuristicOut,", \"SHeuristicOut.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i,i,i);
fclose(HeuristicOut);
                                  


       return 0;
       }
