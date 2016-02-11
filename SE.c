

/* This program generates the model that goes to CPLEX with consideration of spectrum efficiency*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include "SpectrumEfficiency.h"
#include "PowerSaving.h"

#define map_x 50000 //x-axis of map
#define map_y 30000	//y-axis of map
#define Ntp 5100	//number of TPs
#define space 20	//minimum unit on map, should be 20m
#define Rbs 5000	//radius of BS, should be 5000km
#define MNbs 32	//maximum number of possible BS
#define MBR 774	//maximum bandwidth requirement QPSK 3/4

//table:parameters per channel bandwidth
#define BW 10000000	//10 MHz channel bandwidth
#define DSt 5040*1000	//total data subcarriers for DL subframes within 1ms (7 DL symbols per subframe)
// #define DSt 50400*200	//total data subcarriers for 10*200 DL subframe within 5ms*200 (70*200 symbols)
#define PS 120	//pilot subcarriers
#define FFT 1024	//FFT size
#define N 1.12	//n in thermal noise equation
//inputs
#define NF 6	//Rx Noise Figure
#define BH 30	//BS antenna height
#define MH 1	//MS antenna height
#define Gta 16	//BS antenna gain
#define Gto 9	//other DL Tx gain
#define Gra 2	//Rx antenna gain
#define Gro 0	//Rx other gain
#define FFM 3	//fast fading margin

#define MP 20	//UNIT 1mW=30dBm maximum power of a BS: 20W
#define SBP 8.00	//Power for BS in sleep mode
#define ABP 68.73	//Basic power for BS in active mode
#define BP 0.05	//UNIT(%) BLOCKING PROBABILITY
//double R[_t] = {4,6,8,5,2}; // Fill the rates of TPs here

//area of urban and suburban
#define X1urban 750
#define X2urban 1750
#define	Y1urban 500
#define Y2urban 1000
#define X1suburban 500 
#define X2suburban 2000
#define Y1suburban 250
#define Y2suburban 1250
#define UrbanUser 0.7
#define SubUser 0.2
#define RuralUser 0.1

#define LOSmin 20	//los distance
#define NLOSmin 1000	//nlos disance

int Xbs[MNbs+1]; int Ybs[MNbs+1]; // x- and y-coordinates of BS sites
int Xtp[Ntp+1]; int Ytp[Ntp+1]; // x- and y-coordinates or TPs
double D[MNbs+1][Ntp+1];	//distance from BS to TP
double* ptrD;	//pointer to D
double BR[Ntp+1];	//bandwidth requirement of TPs
double* ptrBR;	//pointer to BR

//table:parameters per modulation scheme
int Modulation[Ntp+1];	//Modulation for TP
double Sensitivity[Ntp+1];	//Sensitivity for TP
double DP[Ntp+1];	//Data Bit per Symbol for TP
double SNR[Ntp+1];	//SNR for modulation and TP
int* ptrModulation;	//pointer to Modulation[]
double* ptrSensitivity;	//pointer to Sensitivity[]
double* ptrDP;	//pointer to DP[]
double* ptrSNR;	//pointer to SNR[]

//table:urban corrections
int LT[Ntp+1];	//location type of TP
int* ptrLT;	//pointer to LT[]
int BPL[Ntp+1];	//building penetration loss to TP
int* ptrBPL;	//pointer to BPL

int DS[Ntp+1];	//number of data subcarriers by TP
int* ptrDS;	//pointer to DS[]

double PL[MNbs+1][Ntp+1];	//pathloss between BS and TP
int PathLossType[MNbs+1][Ntp+1];	//pathloss type of TP
double* ptrPL;	//pointer to PL[][]
int* ptrPLT;	//pointer to PLT[]

double P[MNbs+1][Ntp+1];	//power for the transmission from BS to MS
double* ptrP;
 //double Mbt[_t];
 //double Mrr[_r][_r];
 //double Mrt[_r][_t];

void specifyTypes(int Tbs,FILE* PowerSavingCPLEX);
void bounds(int Tbs,FILE* PowerSavingCPLEX);
void checkMSsites(int Tbs);
void FixedCellSize(int Tbs, int* DS, double* ptrD, double* P, double* ptrBR, int* Xbs, int* Ybs, int* Xtp, int* Ytp);
void power(int Tbs, double TN);
void fillCoordinatesBSs(int* Xbs, int* Ybs, int Gbs_x, int Gbs_y, int Nbs_x1, int Nbs_x2, int Nbs_y1, int Nbs_y2, int length);
void fillCoordinatesTPs(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs);
double distance(double xa, double ya, double xb, double yb);
void fillDistance(int* Xbs, int* Ybs, int* Xtp, int* Ytp, double* ptrD, int Tbs);
void fillMs(double* ptrBR);
double ThermalNoise();
void modulation(double TN, double* ptrBR, int* ptrModulation, double* ptrSensitivity, double* ptrDP, double* ptrSNR);
void LocationType(int* ptrLT, int* Xtp, int* Ytp);
void bpl(int* ptrLT, int* ptrBPL);
void DataSubcarriers(int* ptrDS, double* ptrBR, double* ptrDP);
void PathLoss(double* ptrD, int* ptrPLT, int* ptrLT, int Tbs, double* ptrPL);

void constraint9(int Tbs);
void constraint8(int Tbs);
void constraint7(int Tbs);
void constraint6(int Tbs,FILE* PowerSavingCPLEX);
void constraint5(int Tbs,FILE* PowerSavingCPLEX);
void constraint4(int Tbs,FILE* PowerSavingCPLEX);
void constraint3(int Tbs,FILE* PowerSavingCPLEX);
void constraint2(int Tbs,FILE* PowerSavingCPLEX);
void constraint1(int Tbs, int* ptrDS,FILE* PowerSavinCPLEX);
void objective(int Tbs,FILE* PowerSavingCPLEX);
//-------------------------------------------------------
main() {
 int Gbs_x=map_x/Rbs; // Number of grids on y=1
 int Gbs_y=map_y/Rbs; // Number of grids on x=1
 int Nbs_x1=floor(Gbs_x/2);	//Number og BSs on y=1
 int Nbs_x2=floor(((map_x-Rbs)/Rbs)/2);	//Number og BSs on y=2
 int Nbs_y1=floor(Gbs_y/2);	//Number og BSs on x=1
 int Nbs_y2=floor(((map_y-Rbs)/Rbs)/2);	//Number og BSs on x=2
 int Tbs=(Nbs_x1*Nbs_y1)+(Nbs_x2*Nbs_y2);	//Total BS on map
 int length=Rbs/space;	//real length
 int Px=map_x/space;	//number of points on x-axis
 int Py=map_y/space;	//number of points on y-axis
 double TN;	//Thermal Noise
 int i,j;

 ptrD=&D[0][0];
 ptrBR=&BR[0];
 ptrModulation=&Modulation[0];
 ptrSensitivity=&Sensitivity[0];
 ptrDP=&DP[0];
 ptrSNR=&SNR[0];
 ptrBPL=&BPL[0];
 ptrLT=&LT[0];
 ptrDS=&DS[0];
 ptrPLT=&PathLossType[0][0];
 ptrPL=&PL[0][0];
 //ptrP=&P[0][0];
//printf("\n%d\n",Nbs_x1);printf("\n%d\n",Nbs_x2);printf("\n%d\n",Nbs_y1);printf("\n%d\n",Nbs_y2);printf("\n%d\n",Tbs);

 //for(counter=0;counter<_t;counter++)
	 //totalr+=R[counter];//totalr = 4+6+8+5+2


 //checkInputs();
 fillCoordinatesBSs(Xbs,Ybs,Gbs_x,Gbs_y,Nbs_x1,Nbs_x2,Nbs_y1,Nbs_y2,length); 
 /*for(i=1;i<=Tbs;i++)
	 printf("Xbs[%d],Ybs[%d]=(%d,%d)\n",i,i,Xbs[i],Ybs[i]);
 printf("\n");*/
 fillCoordinatesTPs(Xtp,Ytp,Px,Py,Tbs,Xbs,Ybs);
 /*for(i=1;i<=Ntp;i++)
	 printf("Xtp[%d],Ytp[%d]=(%d,%d)\n",i,i,Xtp[i],Ytp[i]);*/
 fillDistance(Xbs,Ybs,Xtp,Ytp,ptrD,Tbs);
 printf("\n");
/*for(i=1;i<=Tbs;i++)
	for(j=1;j<=Ntp;j++)
		printf("D[%d][%d]=%.2lf\n",i,j,D[i][j]);
printf("\n");	*/
 fillMs(ptrBR);
/* for(i=1;i<=Ntp;i++)
	 printf("BR[%d]:%.2lf Mbps\n",i,BR[i]);
 printf("\n");*/
 TN=ThermalNoise();
 //printf("%.2lf\n",TN);
 //printf("\n");
 modulation(TN,ptrBR,Modulation,Sensitivity,DP,SNR);
 //for(i=1;i<=Ntp;i++)
//	 printf("TP[%d] with modulation %d, DP: %.2lf, SNR: %.2lf, S: %.2lf\n",i,Modulation[i],DP[i],SNR[i],Sensitivity[i]);
 //printf("\n");
 LocationType(ptrLT,Xtp,Ytp);
 bpl(ptrLT,ptrBPL);
 //for(i=1;i<=Ntp;i++)
//	 printf("TP[%d]'s location type is %d, bpl is %d\n",i,LT[i],BPL[i]);
 DataSubcarriers(ptrDS,ptrBR,ptrDP);
 //printf("\n");
 //for(i=1;i<=Ntp;i++)
//	 printf("data subcarriers for TP[%d]:%d\n",i,DS[i]);
// printf("\n");
 
 PathLoss(ptrD,ptrPLT,ptrLT,Tbs,ptrPL);
power(Tbs,TN);
FixedCellSize(Tbs,DS,ptrD,P[0],ptrBR,Xbs,Ybs,Xtp,Ytp);

FILE *PowerSavingCPLEX;
if((PowerSavingCPLEX=fopen("PowerSavingCPLEX","w"))==NULL)
	printf("\nerror!Fail to open file!");
else
	printf("\nOpen PowerSavingCPLEX successfully!\n");
fprintf(PowerSavingCPLEX,"This is the input to CPLEX for power saving model.\n");
 objective(Tbs,PowerSavingCPLEX);
fprintf(PowerSavingCPLEX,"st\n");
 printf("st\n");
 constraint1(Tbs,ptrDS,PowerSavingCPLEX);
 constraint2(Tbs,PowerSavingCPLEX);
 constraint3(Tbs,PowerSavingCPLEX);
 constraint4(Tbs,PowerSavingCPLEX);
 constraint5(Tbs,PowerSavingCPLEX);
 //constraint5(nr,nt,totalr);
 constraint6(Tbs,PowerSavingCPLEX);
 constraint7(Tbs);
 constraint8(Tbs);
 constraint9(Tbs);
fprintf(PowerSavingCPLEX,"bounds\n");
 printf("bounds\n");
 bounds(Tbs,PowerSavingCPLEX);

 specifyTypes(Tbs,PowerSavingCPLEX);
fprintf(PowerSavingCPLEX,"end\n");
 printf("end\n");
fclose(PowerSavingCPLEX);
 checkMSsites(Tbs);
	heuristic(Tbs,P[0],Ntp,MP,DSt,BP,DS,Xbs,Ybs,Xtp,Ytp,ptrD,BR);
	printf("===================================================================================================================================\n");
	Sheuristic(Tbs,P[0],Ntp,MP,DSt,BP,DS,Xbs,Ybs,Xtp,Ytp,ptrD,BR); 

FILE *outfile, *outfile1;
 if ((outfile=fopen("outfile1.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen file successfully!\n");
fprintf(outfile,"BS%dMS%dBP%g\n",Tbs,Ntp,BP);
for(i=1;i<=Tbs;i++)
for(j=1;j<=Ntp;j++)
	fprintf(outfile,"%d. BS[%d]=%d,%d MS[%d]=%d,%d DS=%d power=%g mW BW=%g Mbps\n",(i-1)*Ntp+j,i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],DS[j],P[i][j],BR[j]);
fclose(outfile);

if ((outfile1=fopen("coordinates.txt", "w")) == NULL)
printf("\n\nerror!Fail to open file!");
else
printf("\n\nOpen coordinates.txt successfully!\n");

fprintf(outfile1,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
for(i=1;i<=Tbs;i++){
for(j=1;j<=Ntp;j++)
fprintf(outfile1,"%d. BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BW %g Mbps\n",(i-1)*Ntp+j,i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],DS[j],P[i][j],BR[j]);
}
fprintf(outfile1,"\n\n");
for(i=1;i<=Tbs;i++)
fprintf(outfile1," BS[%d]=%d,%d \n",i,Xbs[i],Ybs[i]);

fclose(outfile1);
 //system("pause");
 return;
}
//------------------------------------------------------
//*******************************************************
void checkMSsites(int Tbs){
int i,j;
for(i=1;i<=Tbs;i++)
	for(j=1;j<=Ntp;j++)
		if(i!=j)
			if(Xtp[i]==Xtp[j]&&Ytp[i]==Ytp[j]){
				printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		exit(1);
                }
}
//------------------------------------------------------
//*******************************************************
void power(int Tbs, double TN){
	int i,j;

	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			P[i][j]=pow(10,((-Gta-Gto-Gra-Gro+NF+FFM+BPL[j]+PL[i][j]+Sensitivity[j])/10))/1000;//+TN
}
void PathLoss(double* ptrD, int* ptrPLT, int* ptrLT, int Tbs, double* ptrPL){
	int i,j;

	for(i=1;i<=Tbs;i++){
		for(j=1;j<=Ntp;j++){
			//printf("BS[%d]MS[%d] distance:%.2lf, location type:%d\n",i,j,D[i][j],LT[j]);
			if(LOSmin<=D[i][j]&&D[i][j]<NLOSmin){
				PathLossType[i][j]=1;	//LOS
				PL[i][j]=42.6+26*log10(D[i][j]/1000)+20*log10(BW/1000000);
				//printf("BS[%d]MS[%d] pathloss:%.2lf, type:%d\n",i,j,PL[i][j],PathLossType[i][j]);
			}
			else if(NLOSmin<=D[i][j]){
				PathLossType[i][j]=2;	//NLOS
				if(LT[j]==1)
					PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(BH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+((44.9-(6.55*log10(BH)))*(log10(D[i][j]/1000)))+3;
				else if(LT[j]==2||LT[j]==3)
					PL[i][j]=46.3+(33.9*log10(BW/1000000))-(13.82*log10(BH))-(((1.1*log10(BW/1000000))-0.7)*MH)+((1.56*log10(BW/1000000))-0.8)+(44.9-(6.55*log10(BH)))*(log10(D[i][j]/1000));
				else{
					printf("MS[%d] is in the unregular type.\n");
					exit(1);}
				//printf("BS[%d]MS[%d] pathloss:%.2lf, type:%d\n",i,j,PL[i][j],PathLossType[i][j]);
			}
			else{
				printf("BS[%d] and MS[%d] are at the same spot.\n",i,j);
				exit(1);}
			//printf("BS[%d]MS[%d] DISTANCE:%.2lf, LT:%d, pathloss:%.2lf, type:%d\n",i,j,D[i][j],LT[j],PL[i][j],PathLossType[i][j]);
		}
	}
	return;
}
void FixedCellSize(int Tbs, int* DS, double* ptrD, double* P, double* ptrBR, int* Xbs, int* Ybs, int* Xtp, int* Ytp){
	int i,j;
	double ShortestDistance=1000000000;
	int ConnectionAssignment[Ntp]={-1};
	double TotalBandwidth=0.0, TotalPower=0.0, NOV=0.0;
	int TotalDS=0;
	double Power[MNbs]={0};
	int DataSub[MNbs]={0};
	
	for(j=0;j<Ntp;j++){
		for(i=0;i<Tbs;i++){
			if(D[i+1][j+1]<ShortestDistance){
				ShortestDistance=D[i+1][j+1];
				ConnectionAssignment[j]=i;
			}
		}
		if(ShortestDistance>Rbs){
			printf("%lf MS%d [%d][%d] locates outside the cells!\n Porgram Stops!!\n",ShortestDistance,j+1,Xtp[j+1],Ytp[j+1]);
			exit(1);
		}
		
		ShortestDistance=1000000000;
	}
	
	FILE *FixedCellSize;
 	if ((FixedCellSize=fopen("FixedCellSize.txt", "w")) == NULL)
		printf("\n\nerror!Fail to open file!");
	else
		printf("\n\nOpen FixedCellSize successfully!\n");
/*
	for(i=0;i<Ntp;i++)
		fprintf(FixedCellSize,"MS%d BS%d",i+1,ConnectionAssignment[i]+1);	
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			fprintf(FixedCellSize,"MS%d BS%d Power:%lf\n",i+1,ConnectionAssignment[i]+1),*(P+(i+1)*(Ntp+1)+j+1);
*/
	fprintf(FixedCellSize,"\n#This file is meant to display the regular system.\n");
	fprintf(FixedCellSize,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	for(i=0;i<Tbs;i++)
		fprintf(FixedCellSize," \"BS[%d]\" %d %d %d\n",i+1,Xbs[i+1],Ybs[i+1],250);
	fprintf(FixedCellSize,"\n\n");
        for(i=0;i<Tbs;i++){
		fprintf(FixedCellSize,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
        	for(j=0;j<Ntp;j++){
			if(ConnectionAssignment[i]<0){
				printf("MS has no coonections! Program Stops!!\n");
				exit(1);
			}
			else if(ConnectionAssignment[j]==i){
				fprintf(FixedCellSize,"BS[%d] %d %d BS[%d] %d %d\n",i+1,Xbs[i+1],Ybs[i+1],i+1,Xbs[i+1],Ybs[i+1]);
				fprintf(FixedCellSize,"BS[%d] %d %d MS[%d] %d %d DS %d power %lf mW BR %lf\n",i+1,Xbs[i+1],Ybs[i+1],j+1,Xtp[j+1],Ytp[j+1],DS[j+1],*(P+(i+1)*(Ntp+1)+j+1),BR[j+1]);
			}
		}
	}
	fprintf(FixedCellSize,"\n\n");
	for(i=0;i<Ntp;i++){
		TotalBandwidth+=BR[i+1];
		TotalDS+=DS[i+1];
		for(j=0;j<Tbs;j++)
			if(ConnectionAssignment[i]==j){
				Power[j]+=*(P+(j+1)*(Ntp+1)+i+1);
				if(Power[j]>MP){
					printf("BS%d has too many power demands! Program Stops!\n",j+1);
					exit(1);
				}
				TotalPower+=*(P+(j+1)*(Ntp+1)+i+1);
				DataSub[j]+=*(DS+i+1);
				NOV+=(BR[i+1]/DS[i+1]);
				if(DataSub[j]>DSt){
					printf("BS%d has too many subcarrier demands. Program Stops!\n",j+1);
					exit(1);
				}
			}
	}
	for(i=0;i<Tbs;i++){
		if(Power[i]>0){
			Power[i]+=ABP;
			TotalPower+=ABP;
		}
		else{		
			Power[i]+=SBP;
			TotalPower+=SBP;
		}
	}
	fprintf(FixedCellSize,"Total Power %lf\n",TotalPower);
        fprintf(FixedCellSize,"Total DSs: %d\n",TotalDS);
	fprintf(FixedCellSize,"Total BR: %lf\n",TotalBandwidth);
        //fprintf("Total served MSs: %d\n",TotalServedMSs);
	fprintf(FixedCellSize,"Objective Value: %lf\n",TotalBandwidth/TotalDS/TotalPower);
	fprintf(FixedCellSize,"New Objective Valuee: %lf\n",NOV/TotalPower);
	fprintf(FixedCellSize,"\n\n");
	
	fprintf(FixedCellSize," plot \"FixedCellSize.txt\" index 0:0 using 2:3:1 notitle with labels, \"FixedCellSize.txt\" index 0:0 using 2:3:4 title \"cell size\" with circles, \"coordinates.txt\" index 0:0 using 6:7 title \"users\" with points");
	for(i=1;i<=Tbs;i++)
		fprintf(FixedCellSize,", \"FixedCellSize.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i,i,i);
	fclose(FixedCellSize);
}
void DataSubcarriers(int* ptrDS, double* ptrBR, double* ptrDP){
	int i;

	for(i=1;i<=Ntp;i++)
		DS[i]=ceil((1/N/3*9/8/4*5*FFT)/DP[i]*BR[i]/10);
	//printf("%.2lf\n",1/N/3*9/8*4/5*FFT);
		//BR[i]/(BW*1/N*DP*(1-1/5)*3/9)*8*FFT;
}
void LocationType(int* ptrLT, int* Xtp, int* Ytp){
	int j;

	for(j=1;j<=Ntp;j++){
		if(X1urban<=Xtp[j]&&Xtp[j]<=X2urban&&Y1urban<=Ytp[j]&&Ytp[j]<=Y2urban)
			LT[j]=1;	//urban
		else if(X1suburban<=Xtp[j]&&Xtp[j]<=X2suburban&&Y1suburban<=Ytp[j]&&Ytp[j]<=Y2suburban)
			LT[j]=2;	//suburban
		else
			LT[j]=3;	//rural
	}
	return;
}
void bpl(int* ptrLT, int* ptrBPL){
	int i;
		
	for(i=1;i<=Ntp;i++){
		if(LT[i]==1)
			BPL[i]=18;	//urban
		else if(LT[i]==2)
			BPL[i]=15;	//suburban
		else
			BPL[i]=12;	//rural
	}
	return;
}
void modulation(double TN, double* ptrBR, int* ptrModulation, double* ptrSensitivity, double* ptrDP, double* ptrSNR){

	int i;

	for(i=1;i<=Ntp;i++){
		if(BR[i]<=0.86){
			SNR[i]=6.4;
			Modulation[i]=1;
			DP[i]=0.5;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=1.72){
			SNR[i]=9.4;
			Modulation[i]=2;
			DP[i]=1;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=2.58){
			SNR[i]=11.2;
			Modulation[i]=3;
			DP[i]=1.5;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=3.44){
			SNR[i]=16.4;
			Modulation[i]=4;
			DP[i]=2;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=5.16){
			SNR[i]=18.2;
			Modulation[i]=5;
			DP[i]=3;
			Sensitivity[i]=TN+NF+SNR[i];}
		else if(BR[i]<=6.88){
			SNR[i]=22.7;
			Modulation[i]=6;
			DP[i]=4;
			Sensitivity[i]=TN+NF+SNR[i];}
		else {
			SNR[i]=24.4;
			Modulation[i]=7;
			DP[i]=4.5;
			Sensitivity[i]=TN+NF+SNR[i];}
	}
}
double ThermalNoise(){
	
return	-174+10*log10(BW*N*(720+120+1)/(double)1024);	//DSt+PS
	
//return tn;
}
void specifyTypes(int Tbs,FILE* PowerSavingCPLEX) {
	int i,j;

 	printf("binaries\n");
	fprintf(PowerSavingCPLEX,"binaries\n");
/*
	//dBM
 	for(i=1;i<=Tbs;i++)
  		for(j=1;j<=Ntp;j++)
   			printf("dBM%d_%d \n",i,j);
	
	for(i=1;i<=Tbs;i++)
		printf("dBA%d\n",i);
	for(i=1;i<=Tbs;i++) 
		printf("dBS%d\n",i);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("X%d_%d\n",i,j);
*/
	for(i=1;i<=Tbs;i++) 
		for(j=1;j<=Ntp;j++)
			printf("a%d_%d\n",i,j);//dBMi_j
 	for(i=1;i<=Tbs;i++) 
		printf("b%d\n",i);		
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("c%d_%d \n",i,j);
	
	for(i=1;i<=Tbs;i++)
  		for(j=1;j<=Ntp;j++)
   			fprintf(PowerSavingCPLEX,"dBM%d_%d\n",i,j);	
	for(i=1;i<=Tbs;i++)
		fprintf(PowerSavingCPLEX,"dBA%d\n",i);
/*	for(i=1;i<=Tbs;i++) 
		fprintf(PowerSavingCPLEX,"dBS%d\n",i);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d\n",i,j);*/

return;
}
//*******************************************************
void bounds(int Tbs, FILE* PowerSavingCPLEX) {
	int i, j;

	for(i=1;i<=Tbs;i++) 
		for(j=1;j<=Ntp;j++)
			printf("0 <= u%d_%d\n",i,j);//dBMi_j
 	for(i=1;i<=Tbs;i++) 
		printf("0 <= v%d\n",i);
	//for(i=1;i<=Tbs;i++) 
		printf("0 <= t\n");		
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("0 <= w%d_%d \n",i,j);
/*
	for(i=1;i<=Tbs;i++) 
		for(j=1;j<=Ntp;j++)
			printf("0 <= a%d_%d\n",i,j);//dBMi_j
 	for(i=1;i<=Tbs;i++) 
		printf("0 <= b%d\n",i);		
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("0 <= c%d_%d \n",i,j);
*/
	fprintf(PowerSavingCPLEX,"dBS=1\n");
	//for(i=1;i<=Tbs;i++) 
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= dBM%d_%d\n",i,j);//dBMi_j
 	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBA%d\n",i);
	//for(i=1;i<=Tbs;i++) 
	//	fprintf(PowerSavingCPLEX,"0 <= dBS%d <= 1\n",i);		
	//for(i=1;i<=Tbs;i++)
	//	for(j=1;j<=Ntp;j++)
	//		fprintf(PowerSavingCPLEX,"0 <= X%d_%d\n",i,j);
return;
}
   /*
   //*******************************************************
void checkInputs(void) {
 // In this code, _r should be a perfect square
 double x = sqrt(_r);        // In this code, _r should be a perfect square
 assert(x == floor(x));

 // With _r RS site (let n=sqrt(_r)), there are m=(n-1)^2
 // squares and this m is the limit on TPs. Each square can have one
 // TP in the middle.
 x = sqrt(_r);
 x = pow( (x-1), 2);
 assert(_t <= x);

 // Make sure the arrays TP[] and R[] (of the rates) have a size equal to _t.
 //size_t size_of_array = sizeof(TP) / sizeof(TP[0]);
 //assert(size_of_array == _t);
 //We won't use it since TP is declared as: TP[_t]. It will be true. But hopefully,
 //the user won't forget to fill all _t positions in arrays TP[] and R[].
}
//*******************************************************
*/
void fillDistance(int* Xbs, int* Ybs, int* Xtp, int* Ytp, double* ptrD, int Tbs){
int i,j;
//double k=2;
for(i=1;i<=Tbs;i++)
	for(j=1;j<=Ntp;j++){
		D[i][j]=space*distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j]);
		//printf("%d %d %d %d D[%d][%d]=%.2lf\n",Xbs[i],Ybs[i],Xtp[j],Ytp[j],i,j,D[i][j]);
	}//printf("%.2lf",D[1][1]);
return;
}

void fillCoordinatesTPs(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs){
int i,j;

	double ShortestDistance;
	double tempDistance[Tbs];
	
	for(i=0;i<Tbs;i++)
		tempDistance[i]=0.0;
/*initialize random seed:*/
srand(time(NULL));

/*generate number:*/
for(i=1;i<=ceil(Ntp*RuralUser);i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=1;j<=Tbs;j++){
		tempDistance[j-1]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j-1]<ShortestDistance)
			ShortestDistance=tempDistance[j-1];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		printf("%d renew\n",i);}
}

for(i=ceil(Ntp*RuralUser);i<=ceil(Ntp*(RuralUser+SubUser));i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2suburban-X1suburban+1)+X1suburban;
	Ytp[i]=rand()%(Y2suburban-Y1suburban+1)+Y1suburban;
	for(j=1;j<=Tbs;j++){
		tempDistance[j-1]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j-1]<ShortestDistance)
			ShortestDistance=tempDistance[j-1];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		printf("%d renew\n",i);}
}
for(i=ceil(Ntp*(RuralUser+SubUser));i<=Ntp;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(X2urban-X1urban+1)+X1urban;
	Ytp[i]=rand()%(Y2urban-Y1urban+1)+Y1urban;
	for(j=1;j<=Tbs;j++){
		tempDistance[j-1]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j-1]<ShortestDistance)
			ShortestDistance=tempDistance[j-1];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		printf("%d renew\n",i);}
}/*
for(i=41;i<=64;i++){
	Xtp[i]=(rand()%(376))+250;
	Ytp[i]=(rand()%(376))+250;
}*/
return;
}

void fillCoordinatesBSs(int* Xbs, int* Ybs, int Gbs_x, int Gbs_y, int Nbs_x1, int Nbs_x2, int Nbs_y1, int Nbs_y2, int length){
int Sx=1;
int Sy=1;
int Zy=Nbs_x1;
int Ay=Nbs_x1;
int By=Nbs_x2;
int Zx=Nbs_y1;
int Ax=Nbs_y1;
int Bx=Nbs_y2;
int i,j,k,l;

 //n = sqrt(_r); // number of RS sites along an axis of the grid.
 //printf("n:%.0lf\n",n)
 
 for(j=1; j<=(Gbs_y-1); j++) {
	 for(i=Sy; i<=Zy; i++)
		 Ybs[i] = j*length;
	 Sy=Zy+1;//printf("%d,%d\n",Sy,Zy);
	 if(j%2==1)
	 Zy+=By;
	 else
	 Zy+=Ay;	 
 }
 //printf("%d,%d",Xbs[4],Ybs[4]);
 //printf("%d",Tbs);
  //printf("Xrs[%d]=%lf, Yrs[%d]=%lf\n",i,Xrs[i],i,Yrs[i]);  
 for(k=1; k<=(Gbs_x-1); k++){
	for(l=Sx; l<=Sx+(Ay+By)*(Zx-1); l=l+Ay+By)
		Xbs[l]=k*length;
	if(k%2==1){
		Sx+=Ay;
		Zx+=Bx;
	}
	else{
		Sx=Sx-Ay+1;
		Zx+=Ax;
	}
 }
// printf("%d,%d",Xbs[4],Ybs[4]);
 
 //Xbs = (double)(n-1)/2;
 //Ybs = 0;
 //printf("Xbs=%lf, Ybs=%lf\n",Xbs,Ybs);

 //for(i=0; i<_t; i++) {
  //Xtp[i] = ((int)TP[i]%(int)(n-1))+0.5;
  //Ytp[i] = floor( TP[i]/(n-1) )+0.5;
  //printf("Xtp[%d]=%lf, Ytp[%d]=%lf\n",i,Xtp[i],i,Ytp[i]);
 //}
 return;
}
//*******************************************************

double distance(double xa, double ya, double xb, double yb) {
 double value;
 //double temp1, temp2;

// temp1 = xa-xb;
// temp1 = pow(temp1,2);
// temp2 = ya-yb;
// temp2 = pow(temp2,2);
// value = sqrt(temp1+temp2);

 value = sqrt( pow((xa-xb),2) + pow((ya-yb),2) );

 return value;
}
//*******************************************************
void fillMs(double* ptrBR) { // Mbr,Mbt,Mrr,Mrt/* The grid is a square with length=width. The BS is at the top line in the middle.
/*length is even, the BS will coincide with a RS site. Then we "cancel" this RS,
which means mbr=0. See function_fillMbr.ps illustration.*/
	int i;
	//double dist;
	double j;
	srand(time(NULL));

	for(i=1;i<=Ntp;i++){
		//BR[i]
		j=(rand()%(MBR+1));
		while(j==0)
			j=(rand()%(MBR+1));
		BR[i]=j/100;
		//printf("%.2lf\n",j);
		//printf("%.2lf\n",BR[i]);
	}
//	for(i=0;i<_r;i++) {
  //dist = distance(Xbs,Ybs,Xrs[i],Yrs[i]); //BS colocated with RS 
  //if(dist==0) Mbr[i]=0; 
  //else {  
	//  dist = ceil(dist);  
  //if(dist==1) 
	//  Mbr[i] = 10; //10 Mbps   
  //else if(dist==2) Mbr[i] = 5; // 5 Mbps
   //else if(dist==3) Mbr[i] = 2;   
   //else if(dist==4) Mbr[i] = 1;  
   //else Mbr[i] = 0;
  //} //printf("Mbr[%d]=%lf, ",i,Mbr[i]);
  //if(i%(int)sqrt(_r)==sqrt(_r)-1) printf("\n");
 //}
//	for(i=0;i<_t;i++) {
//		dist = distance(Xbs,Ybs,Xtp[i],Ytp[i]);
  //dist = ceil(dist); 
  //if(dist==1) Mbt[i]=10;
  //else if(dist==2) Mbt[i]=5;
  //else if(dist==3) Mbt[i]=2;
  //else if(dist==4) Mbt[i]=1;
  //else Mbt[i]=0;//printf("Mbt of TP %.0lf = %lf\n", TP[i], Mbt[i])
 //}
	
//	for(i=0;i<_r;i++) {
  //for(j=0;j<_r;j++){
   //dist = distance(Xrs[i],Yrs[i],Xrs[j],Yrs[j]);
   //if(dist==0) Mrr[i][j]=0;   else{   dist = ceil(dist);
    //if(dist==1) Mrr[i][j]=10;
    //else if(dist==2) Mrr[i][j]=5;
    //else if(dist==3) Mrr[i][j]=2;
    //else if(dist==4) Mrr[i][j]=1;
    //else Mrr[i][j]=0;
   //}   //printf("RS%d - RS%d, distance=%lf, Mrr=%lf\n",i,j,dist,Mrr[i][j]);
  //}
  //printf("\n");} for(i=0;i<_r;i++) {
  //for(j=0;j<_t;j++){   dist = distance(Xrs[i],Yrs[i],Xtp[j],Ytp[j]);
   //assert(dist != 0); // Since TP are inside squared229,1         51%dist = distance(Xrs[i],Yrs[i],Xtp[j],Ytp[j]);
   //assert(dist != 0); // Since TP are inside squared
   //dist = ceil(dist);
   //if(dist==1) Mrt[i][j]=10;
   //else if(dist==2) Mrt[i][j]=5;
   //else if(dist==3) Mrt[i][j]=2;
   //else if(dist==4) Mrt[i][j]=1;
   //else Mrt[i][j]=0;
   //printf("RS%d - TP%d, distance=%lf, Mrt = %lf\n",i,j,dist,Mrt[i][j]);
  //}
  //printf("\n");
 //}

 return;
}
//*******************************************************
void constraint9(int Tbs) {
 int i,j;

	for(i=1;i<=Tbs;i++){
		printf("%lf v%d + ",ABP-SBP,i);
		for(j=1;j<=Ntp;j++){
			if(i==Tbs&&j==Ntp)
				printf("%lf u%d_%d",P[i][j],i,j);
			else
				printf("%lf u%d_%d + ",P[i][j],i,j);
		}
	}
	printf(" + %lf t = 1\n",Tbs*SBP);
 /*for(i=0; i<=nr; i++)
 printf("fbr%d <= %.2lf\n",i,Mbr[i]);

 for(j=0; j<=nt; j++)
 printf("fbt%d <= %.2lf\n", j,Mbt[j]);

 for(i=0; i<=nr; i++)
 for(j=0; j<=nr; j++)
 printf("frr%d_%d <= %.2lf\n", i,j,Mrr[i][j]);

 for(i=0; i<=nr; i++)
 for(j=0; j<=nt; j++)
 printf("frt%d_%d <= %.2lf\n", i,j,Mrt[i][j]);*/

 return;
}
//*******************************************************
void constraint8(int Tbs) {
 int i,j;
	int Q2=1002;

	for(i=1;i<=Tbs;i++){
		printf("v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		printf("v%d - t - %d b%d <= %d\n",i,Q2,i,Q2);
		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d >= -%d\n",i,Q2,i,Q2);
		//fprintf(PowerSavingCPLEX,"v%d - t - %d b%d <= %d\n",i,Q2,i,Q2);
		for(j=1;j<=Ntp;j++){
			printf("u%d_%d - t - %d a%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
			printf("u%d_%d - t - %d a%d_%d <= %d\n",i,j,Q2,i,j,Q2);
			printf("w%d_%d - t - %d c%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
			printf("w%d_%d - t - %d c%d_%d <= %d\n",i,j,Q2,i,j,Q2);
			//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
			//fprintf(PowerSavingCPLEX,"u%d_%d - t - %d a%d_%d <= %d\n",i,j,Q2,i,j,Q2);
			//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d >= -%d\n",i,j,Q2,i,j,Q2);
			//fprintf(PowerSavingCPLEX,"w%d_%d - t - %d c%d_%d <= %d\n",i,j,Q2,i,j,Q2);
		}
	}
		
 /*for(i=0; i<=nt; i++) {
  printf("y%d",i);

  for(j=0; j<=nr; j++)
  printf("+w%d_%d",j,i);

  printf("=%.2lf\n",R[i]);
 }*/
 return;
}
//*******************************************************
void constraint7(int Tbs) {
	int i,j;
	int Q1=1001;

/*	for(i=1; i<=Tbs; i++)
		printf("dBA%d + dBS%d = 1\n",i,i);
	for(i=1; i<=Tbs; i++)
		fprintf(PowerSavingCPLEX,"dBA%d + dBS%d = 1\n",i,i);*/
	for(i=1;i<=Tbs;i++){
		printf("v%d - %d b%d <= 0\n",i,Q1,i);
		//fprintf(PowerSavingCPLEX,"v%d <= %d t\n",i,Q1);
		for(j=1;j<=Ntp;j++){
			printf("u%d_%d - %d a%d_%d <= 0\n",i,j,Q1,i,j);
			printf("w%d_%d - %d c%d_%d <= 0\n",i,j,Q1,i,j);
			//fprintf(PowerSavingCPLEX,"u%d_%d <= %d t\n",i,j,Q1);
			//fprintf(PowerSavingCPLEX,"w%d_%d <= %d t\n",i,j,Q1);
		}
	}
return;
}
//*******************************************************
void constraint6(int Tbs,FILE* PowerSavingCPLEX) {
	int i,j;
	int Q=1000;

	for(i=1;i<=Tbs;i++){
		for(j=1;j<Ntp;j++)
			printf("w%d_%d - u%d_%d + ",i,j,i,j);
		printf("w%d_%d - u%d_%d = 0\n",i,j,i,j);
	}
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("w%d_%d - %d v%d - u%d_%d + %d t >= 0\n",i,j,Q,i,i,j,Q);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("w%d_%d - u%d_%d <= 0\n",i,j,i,j);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("w%d_%d >= 0\n",i,j);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			printf("w%d_%d - %d v%d <= 0\n",i,j,Q,i);
	for(i=1;i<=Tbs;i++){
		for(j=1;j<Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d - dBM%d_%d + ",i,j,i,j);
		fprintf(PowerSavingCPLEX,"X%d_%d - dBM%d_%d = 0\n",i,j,i,j);
	}
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d - %d dBA%d - dBM%d_%d >= -%d\n",i,j,Q,i,i,j,Q);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d - dBM%d_%d <= 0\n",i,j,i,j);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d >= 0\n",i,j);
	for(i=1;i<=Tbs;i++)
		for(j=1;j<=Ntp;j++)
			fprintf(PowerSavingCPLEX,"X%d_%d - %d dBA%d <= 0\n",i,j,Q,i);
 return;
}
//*******************************************************
void constraint5(int Tbs,FILE* PowerSavingCPLEX) {
	int i,j;

 	for(i=1;i<=Tbs;i++)
	 	for(j=1;j<=Ntp;j++)
			if(j==Ntp)
				printf("u%d_%d - v%d >= 0\n",i,j,i);
 			else
				printf("u%d_%d + ",i,j);
	for(i=1;i<=Tbs;i++)
	 	for(j=1;j<=Ntp;j++)
			if(j==Ntp)
				fprintf(PowerSavingCPLEX,"dBM%d_%d - dBA%d >= 0\n",i,j,i);
 			else
				fprintf(PowerSavingCPLEX,"dBM%d_%d + ",i,j);
return;
}
//*******************************************************
void constraint4(int Tbs,FILE* PowerSavingCPLEX) {
 int i,j;

 for(i=1; i<=Tbs; i++) {
  //printf("dbt%d",i);
	 for(j=1; j<=Ntp; j++){
	  if(j==Ntp&&i==Tbs)
		  printf("u%d_%d",i,j);
	else 
		printf("u%d_%d +",i,j);
	 }
  
 }
printf(" - %.2lf t >= 0\n",Ntp*(1-BP));
	for(i=1; i<=Tbs; i++) {
  //printf("dbt%d",i);
	 for(j=1; j<=Ntp; j++){
	  if(j==Ntp&&i==Tbs)
		  fprintf(PowerSavingCPLEX,"dBM%d_%d",i,j);
	else 
		fprintf(PowerSavingCPLEX,"dBM%d_%d +",i,j);
	 }
  
 }
fprintf(PowerSavingCPLEX," >= %.2lf\n",Ntp*(1-BP));
 return;
}
//*******************************************************
void constraint3(int Tbs,FILE* PowerSavingCPLEX) {
 int i,j;
	 for(i=1; i<=Ntp; i++)
		 for(j=1; j<=Tbs; j++) {
	  if(j==Tbs)
		  printf("u%d_%d - t <= 0\n",j,i);
	  else
		  printf("u%d_%d +",j,i);}
		 //printf("<=1\n");
	for(i=1; i<=Ntp; i++)
		 for(j=1; j<=Tbs; j++) {
	  if(j==Tbs)
		  fprintf(PowerSavingCPLEX,"dBM%d_%d <= 1\n",j,i);
	  else
		  fprintf(PowerSavingCPLEX,"dBM%d_%d +",j,i);}
 return;
}

//*******************************************************
void constraint2(int Tbs,FILE* PowerSavingCPLEX) {
	int i,j;

		for(i=1; i<=Tbs; i++){
		for(j=1;j<=Ntp;j++){
			if(j==Ntp)
				printf("%lf u%d_%d", P[i][j],i, j);
			else 
				printf("%lf u%d_%d + ",P[i][j], i, j);
		}
	printf(" - %d t <= 0\n",MP);
	}
	for(i=1; i<=Tbs; i++){
		for(j=1;j<=Ntp;j++){
			if(j==Ntp)
				fprintf(PowerSavingCPLEX,"%lf dBM%d_%d", P[i][j],i, j);
			else 
				fprintf(PowerSavingCPLEX,"%lf dBM%d_%d + ",P[i][j], i, j);
		}
	fprintf(PowerSavingCPLEX," <= %d\n",MP);
	}
 return;
}//*******************************************************

void constraint1(int Tbs, int* ptrDS, FILE* PowerSavingCPLEX) {
	int i,j;
	for(i=1; i<=Tbs; i++){
		for(j=1;j<=Ntp;j++){
			if(j==Ntp)
				printf("%d u%d_%d", DS[j],i, j);
			else 
				printf("%d u%d_%d + ",DS[j], i, j);
		}
	printf(" - %d t <= 0\n",DSt);	//DSt*200 for 1 second
	}
	for(i=1; i<=Tbs; i++){
		for(j=1;j<=Ntp;j++){
			if(j==Ntp)
				fprintf(PowerSavingCPLEX,"%d dBM%d_%d", DS[j],i, j);
			else 
				fprintf(PowerSavingCPLEX,"%d dBM%d_%d + ",DS[j], i, j);
		}
	fprintf(PowerSavingCPLEX," <= %d\n",DSt);	//DSt*200 for 1 second
	}
 return;
}//*******************************************************

void objective(int Tbs, FILE* PowerSavingCPLEX) {
	int count,i;
 
	printf("maximize\n");
	fprintf(PowerSavingCPLEX,"minimize\n");
	for(count=1;count<=Ntp;count++)
		for(i=1;i<=Tbs;i++){
			if(count==Ntp&&i==Tbs)
				printf("%lf u%d_%d\n",BR[count]/DS[count],i,count);
			else 
				printf("%lf u%d_%d + ",BR[count]/DS[count],i,count);
		}
	/*
	for(count=1; count<=Ntp; count++){
		if(count==1)
			printf("((%lf ",BR[count]);
		else
			printf("(%lf ",BR[count]);
		for(i=1;i<=Tbs;i++){
			if(i==1)
				printf("(dBM%d_%d + ",i,count);
			else if(i==Tbs&&count==Ntp)
				printf("dBM%d_%d)))",i,count);
			else if(i==Tbs)
				printf("dBM%d_%d)) + ",i,count);
			else
				printf("dBM%d_%d + ",i,count);
		}
	}
	*/
	/*printf("-");
	for(count=1;count<=Ntp;count++)
		for(i=1;i<=Tbs;i++){
			if(count==Ntp&&i==Tbs)
				printf("%d dBM%d_%d",DS[count],i,count);
			else 
				printf("%d dBM%d_%d - ",DS[count],i,count);
		}*/
	for(count=1; count<=Tbs; count++)
		for(i=1;i<=Ntp;i++){
			if(count==1&&i==1)
				fprintf(PowerSavingCPLEX,"%lf dBS + %lf dBM%d_%d + ",SBP*Tbs,P[count][i],count,i);
			else if(count!=Tbs&&i==Ntp){
				fprintf(PowerSavingCPLEX,"%lf dBM%d_%d",P[count][i],count,i);
				fprintf(PowerSavingCPLEX," + %lf dBA%d + ",ABP-SBP,count);
				//fprintf(PowerSavingCPLEX," + %lf + ",SBP);
			}
			else if(count==Tbs&&i==Ntp){
				fprintf(PowerSavingCPLEX,"%lf dBM%d_%d",P[count][i],count,i);
				fprintf(PowerSavingCPLEX," + %lf dBA%d\n",ABP-SBP,count);
				//fprintf(PowerSavingCPLEX," + %lf\n",SBP*Tbs);
			} 
			else
				fprintf(PowerSavingCPLEX,"%lf dBM%d_%d + ",P[count][i],count,i);
		}
	/*printf("-");
 	for(count=1; count<=Tbs; count++)
		for(i=1;i<=Ntp;i++){
			if(count==1&&i==1)
				printf("%lf dBM%d_%d - ",P[count][i],count,i);
			else if(count!=Tbs&&i==Ntp){
				printf("%lf dBM%d_%d",P[count][i],count,i);
				printf(" - %lf dBA%d",ABP,count);
				printf(" - %lf dBS%d  - ",SBP,count);
			}
			else if(count==Tbs&&i==Ntp){
				printf("%lf dBM%d_%d",P[count][i],count,i);
				printf(" - %lf dBA%d",ABP,count);
				printf(" - %lf dBS%d \n",SBP,count);
			} 
			else
				printf("%lf dBM%d_%d - ",P[count][i],count,i);
		}*/
return;
}
