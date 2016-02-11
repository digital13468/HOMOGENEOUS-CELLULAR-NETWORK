#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"CoverageComputation.h"
#include<math.h>

#define Tbs 11
#define Ntp 5100
#define BP 0.05
#define SBP 8.00	//Power for BS in sleep mode
#define ABP 68.73	//Basic power for BS in active mode

int Xb[Tbs]={0};
int Yb[Tbs]={0};
int B[Tbs]={0};
int Xtp[Ntp]={0};
int Ytp[Ntp]={0};
int dBM[Tbs][Ntp];
int dBA[Tbs]={0};
int A[Tbs][Ntp];
int C[Tbs][Ntp];
double P[Tbs][Ntp];
double BR[Ntp];
double u[Tbs][Ntp];
double w[Tbs][Ntp];
int DataSub[Ntp];
double CumulativeP[Tbs]={0};
double CumulativeBR[Tbs]={0};
double v[Tbs]={0};
int CumulativeDS[Tbs]={0};
int CumSerMS[Tbs]={0};
double TotalP=0.0;
double TotalBR=0.0;
double t=0.0;
double NOV=0.0;
int TotalDS=0;
int TotSerMS=0;

InfoPrinter(FILE* CplexOut){
               int i,j;
               for(i=0;i<Tbs;i++)
               fprintf(CplexOut,"       BS%3d",i+1);
               fprintf(CplexOut,"\n");
               
               for(i=0;i<Ntp;i++){
                                  fprintf(CplexOut,"MS%4d%4d",i+1,dBM[0][i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"      %6d",dBM[j][i]);
                                  fprintf(CplexOut,"\n");
                                  
                                  fprintf(CplexOut," Power%10.5f",P[0][i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"  %10.5f",P[j][i]);
                                  fprintf(CplexOut,"\n");
                                  
                                  fprintf(CplexOut," DS %6d",DataSub[i]);
                                  for(j=1;j<Tbs;j++)
                                  fprintf(CplexOut,"      %6d",DataSub[i]);
                                  fprintf(CplexOut,"\n");
                                  }
}

double * ConnectionPrinter(FILE* CplexOut){
               int i,j;
	double *coverage;
	
	coverage=CoverageComputation(dBM[0],Tbs,Ntp,Xb,Yb,Xtp,Ytp);

fprintf(CplexOut,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);

fprintf(CplexOut,"# plot \"outfile.txt\" index 0:0 using 2:3:1 notitle with labels, \"outfile.txt\" index 0:0 using 2:3:4 title \"cell size\" with circles, \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
for(i=1;i<=Tbs;i++)
fprintf(CplexOut,", \"outfile.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i,i,i);
fprintf(CplexOut,"\n\n");

for(i=0;i<Tbs;i++)
fprintf(CplexOut," \"BS[%d]\" %d %d %g\n",i+1,Xb[i],Yb[i],*(coverage+i));
fprintf(CplexOut,"\n\n");
               for(i=0;i<Tbs;i++){
fprintf(CplexOut,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xb[i],Yb[i],i+1,Xb[i],Yb[i]);
               for(j=0;j<Ntp;j++)
		if(dBM[i][j]==1){
fprintf(CplexOut,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xb[i],Yb[i],i+1,Xb[i],Yb[i]);
fprintf(CplexOut,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW SE %g\n",i+1,Xb[i],Yb[i],j+1,Xtp[j],Ytp[j],DataSub[j],P[i][j],BR[j]/DataSub[j]/P[i][j]);
}

fprintf(CplexOut,"\n\n");
}
}
                                  
main(){
       
FILE *infile, *outfile; 
int i,j;
char id[100];
char tempC;
int a, b, c;
int BS,MS,value;
double fvalue;
size_t loc;
char string[12]="<variables>" ;
char string1[13]="</variables>" ;
char string2[9]="<variable";

if ((infile=fopen("SECRD.txt", "r")) == NULL)
	printf("Fail to open file!");
else
	printf("Open file successfully!\n");
/*if((outfile=fopen("outfile.txt","w"))==NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");*/

while (fscanf(infile, "%s", &id)!=EOF){
	if (strcmp(id,string)==0){
                          fscanf(infile, "%s", &id);
                         while(strcmp(id,string1)!=0)
                         {//printf("Outer Scanning\n");
                                                             while(strcmp(id,string2)==0){
								fscanf(infile, " name=\"");
								fscanf(infile, "%c", &tempC);
								//printf("%c\n",tempC);
								if(tempC=='a'){
									fscanf(infile, "%d_%d\"", &BS, &MS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%d\"/>", &value);
									A[BS-1][MS-1]=value;
									//fscanf(infile, "%s",&id);
									//printf("%s\n",id);
									//fscanf(infile, "%s",&id);
									//printf("%s\n",id);
									//fscanf(infile, "%s",&id);
									//printf("%s\n",id);
								}
								else if(tempC=='b'){
									fscanf(infile, "%d\"", &BS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%d\"/>", &value);
									B[BS-1]=value;
								}
								else if(tempC=='c'){
									fscanf(infile, "%d_%d\"", &BS, &MS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%d\"/>", &value);
									C[BS-1][MS-1]=value;
									
								}
								else if(tempC=='u'){
									fscanf(infile, "%d_%d\"", &BS, &MS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%lf\"/>", &fvalue);
									//printf("%lf\n",fvalue);
									u[BS-1][MS-1]=fvalue;
									//printf("%lf\n",fvalue);
								}
								else if(tempC=='v'){
									fscanf(infile, "%d\"", &BS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%lf\"/>", &fvalue);
									v[BS-1]=fvalue;
								}
								else if(tempC=='w'){
									fscanf(infile, "%d_%d\"", &BS, &MS);
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%lf\"/>", &fvalue);
									w[BS-1][MS-1]=fvalue;
								}
								else if(tempC=='t'){
									fscanf(infile, "\"");
									fscanf(infile, "%s", &id);
									fscanf(infile, " value=\"%lf\"/>", &fvalue);
									t=fvalue;
								}
								/*else if(tempC=='d'){
									fscanf(infile, "%c%c", &tempC,&tempC);
									//printf("%c\n",tempC);
									if(tempC=='M'){
										fscanf(infile, "%d_%d\"", &BS, &MS);
										fscanf(infile, "%s", &id);
                                                             			fscanf(infile, " value=\"%d\"/>", &value);
                                                             			dBM[BS-1][MS-1]=value;
										printf("%d %d\n",BS,MS);
									}
									else{
										fscanf(infile, "%s",&id);
										//printf("%s\n",id);
										fscanf(infile, "%s",&id);
										//printf("%s\n",id);
										fscanf(infile, "%s",&id);
										//printf("%s\n",id);
									}
								}*/
								fscanf(infile, "%s", &id);
                                                             /*fscanf(infile, " name=\"dBM%d_%d\"", &BS, &MS);
                                                             fscanf(infile, "%s", &id);
                                                             fscanf(infile, " value=\"%d\"/>", &value);
                                                             dBM[BS-1][MS-1]=value;
                                                             fscanf(infile, "%s", &id);
								//printf("Inner Scanning\n");
								if(BS==Tbs&&MS==Ntp)
									strcpy(id,string1);*/
                                                             }                                                          
                         }                                                            
	}

}

/*for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
fprintf(outfile,"dBM%d_%d=%d\n",i+1,j+1,dBM[i][j]);*/
fclose(infile);
//fclose(outfile);

FILE *CplexOut;
FILE *CplexIn;
//int BS,MS,value;
int Xbs, Ybs, Xms, Yms, DS;
double p,br;
char id1[100];

if ((CplexIn=fopen("outfile1.txt", "r")) == NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");
if ((CplexOut=fopen("SpectrumEfficiencyCPLEX_Output.txt", "w")) == NULL)
printf("Fail to open file!");
else
printf("Open file successfully!\n");
fscanf(CplexIn, "%s", &id1);
while (fscanf(CplexIn, "%s", &id1)!=EOF){
      //fprintf(CplexOut, "%s\n",id1);
      fscanf(CplexIn, " BS[%d]=%d,%d", &BS, &Xbs, &Ybs);
	Xb[BS-1]=Xbs;
	Yb[BS-1]=Ybs;
      //fprintf(CplexOut,"BS[%d]=(%d,%d) ", BS, Xbs, Ybs);
      fscanf(CplexIn, " MS[%d]=%d,%d", &MS, &Xms, &Yms);
	Xtp[MS-1]=Xms;
	Ytp[MS-1]=Yms;
      //fprintf(CplexOut,"MS[%d]=(%d,%d) ", MS, Xms, Yms);
      fscanf(CplexIn, " DS=%d", &DS);
      DataSub[MS-1]=DS;
      //fprintf(CplexOut,"DS:%d ", DS);
      fscanf(CplexIn, " power=%lf mW", &p);
      P[BS-1][MS-1]=p;
	fscanf(CplexIn," BW=%lf",&br);
	BR[MS-1]=br;
      //fprintf(CplexOut,"power:%g\n", p);
      fscanf(CplexIn," %s",&id1);

}
/*for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
fprintf(outfile,"dBM%d_%d=%d P:%g DS:%d\n",i+1,j+1,dBM[i][j],P[i][j],DataSub[j]);*/
fclose(CplexIn);

for(i=0;i<Tbs;i++){
	dBA[i]=v[i]/t;
	for(j=0;j<Ntp;j++){
		//fprintf(CplexOut,"%d %d %lf %lf %lf\n",i,j,u[i][j], t, u[i][j]/t);
		dBM[i][j]=u[i][j]/t;
	}
}
		
for(i=0;i<Tbs;i++){
	if(dBA[i]!=B[i]){
		printf("dBA[%d]=%d, B[%d]=%d. Program Stops!\n",i+1,dBA[i],i+1,B[i]);
		exit(1);
	}
	for(j=0;j<Ntp;j++){
		if(dBM[i][j]!=A[i][j]){
			printf("dBM[%d][%d]=%d, A[%d][$d]=%d. Program Stops!\n",i+1,j+1,dBM[i][j],i+1,j+1,B[i]);
			exit(1);
		}
	}	
}
for(i=0;i<Tbs;i++)
for(j=0;j<Ntp;j++)
if(dBM[i][j]==1){
CumulativeP[i]+=P[i][j];
TotalP+=P[i][j];
CumulativeDS[i]+=DataSub[j];
	CumulativeBR[i]+=BR[j];
TotalDS+=DataSub[j];
	TotalBR+=BR[j];
CumSerMS[i]++;
TotSerMS++;
	NOV+=(BR[j]/DataSub[j]);
}
for(i=0;i<Tbs;i++){
	if(CumulativeP[i]>0){
		CumulativeP[i]+=ABP;
		TotalP+=ABP;
	}
	else{
		CumulativeP[i]+=SBP;
		TotalP+=SBP;
	}
}

ConnectionPrinter(CplexOut);
InfoPrinter(CplexOut);

for(i=0;i<Tbs;i++){
	//printf("Power from BS %d: %lf\n",i+1,CumulativeP[i]);
       fprintf(CplexOut,"Power from BS %d: %lf\n",i+1,CumulativeP[i]);
       fprintf(CplexOut,"DS from BS %d: %d\n",i+1,CumulativeDS[i]);
	fprintf(CplexOut,"Bandwidth from BS %d: %lf\n",i+1,CumulativeBR[i]);
       fprintf(CplexOut,"MSs served by BS %d: %d\n\n",i+1,CumSerMS[i]); 
       }
       fprintf(CplexOut,"Total Power %lf\n",TotalP);
       fprintf(CplexOut,"Total DSs: %d\n",TotalDS);
       fprintf(CplexOut,"Total served MSs: %d\n",TotSerMS);
	fprintf(CplexOut,"Total Bandwidth: %lf\n",TotalBR);
	fprintf(CplexOut,"Objective Value: %lf\n",TotalBR/TotalDS/TotalP);
	fprintf(CplexOut,"New Objective Value: %lf\n",NOV/TotalP);
fclose(CplexOut);
//system("pause");
}

