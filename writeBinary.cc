#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "CalTV.h"
#include "Wepl.h"
//#include <TProject3D.h>
//ccompile by: g++ -c `root-config --cflags` Tracker.C  ; g++ -o Tracker.out Tracker.o `root-config --glibs`
using namespace std;

double str8Line(double t, double par[])
   {
      double f = par[0] +t*par[1];
      return f;
   }

double str8LineInv(double z, double par1[])
   {
      double fInv = (z-par1[0])/par1[1];
      return fInv;
   }

/////////////////////////////////////////Straight line to calculate T,V coordinates for the TV covoid
void GetTVcoord (double Zc, double &Tc,double &Vc, float plane3[], float plane4[])
{
//	float Z=216.9 + 40; // approx position for the calorimeter entrance
	double ParamZ[2]={plane3[2], plane4[2]-plane3[2]};
	double ParamX[2]={plane3[0], plane4[0]-plane3[0]};
	double ParamY[2]={plane3[1], plane4[1]-plane3[1]};

	double T = str8LineInv(Zc,ParamZ);
	double X = str8Line(T,ParamX);
	double Y = str8Line(T,ParamY);
	Tc=X;
	Vc=Y;
}
/////////////////////////////////////////////////////
///////////////////////////////////////////write binary
void WriteBinaryFile(float projection_angle, char DATA_SOURCE[], int event_counter, float Vw0[],float Vw1[],float Vw2[],float Vw3[], float Tw0[],float Tw1[],float Tw2[],float Tw3[], float uw0[],float uw1[],float uw2[],float uw3[], float calculated_WEPL[],char filenm[])
{
//cout<< "WriteHeader "<<WriteHeader<<" projection_angle "<< projection_angle<< endl;
	ofstream data_file;
	data_file.open(filenm,  ios::binary |ios::trunc);
    char magic_number[] = "PCTD";
	char PHANTOM_NAME[]= "Sensitom";
	char PREPARED_BY[]= "Pier";
//	float projection_angle=0;
	float beam_energy= 200.;
    int version_id = 0;
    int current_time = time(NULL);
    int phantom_name_size = sizeof(PHANTOM_NAME);
    int data_source_size = sizeof(DATA_SOURCE);
    int prepared_by_size = sizeof(PREPARED_BY);
    // Write headers:
	    data_file.write(magic_number, 4); // magic number identifier (note that it doesn't include null terminator '\0')
	    data_file.write(reinterpret_cast<char*>(&version_id), sizeof(int)); // format version identifier
	    data_file.write(reinterpret_cast<char*>(&event_counter), sizeof(int)); // number of events in file
	    data_file.write(reinterpret_cast<char*>(&projection_angle), sizeof(float)); // projection angle
	    data_file.write(reinterpret_cast<char*>(&beam_energy), sizeof(float)); // beam energy
	    data_file.write(reinterpret_cast<char*>(&current_time), sizeof(int)); // generation date
	    data_file.write(reinterpret_cast<char*>(&current_time), sizeof(int)); // pre-process date
	    data_file.write(reinterpret_cast<char*>(&phantom_name_size), sizeof(int));
	    data_file.write(PHANTOM_NAME, phantom_name_size); // phantom name or description (string)
	    data_file.write(reinterpret_cast<char*>(&data_source_size), sizeof(int));
	    data_file.write(DATA_SOURCE, data_source_size); // data source (string)
	    data_file.write(reinterpret_cast<char*>(&prepared_by_size), sizeof(int));
	    data_file.write(PREPARED_BY, prepared_by_size); // prepared by (string)

    // Write event data:
    data_file.write(reinterpret_cast<char*>(Tw0), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Tw1), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Tw2), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Tw3), event_counter * sizeof(float));

    data_file.write(reinterpret_cast<char*>(Vw0), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Vw1), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Vw2), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(Vw3), event_counter * sizeof(float));

    data_file.write(reinterpret_cast<char*>(uw0), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(uw1), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(uw2), event_counter * sizeof(float));
    data_file.write(reinterpret_cast<char*>(uw3), event_counter * sizeof(float));

    data_file.write(reinterpret_cast<char*>(calculated_WEPL), event_counter * sizeof(float));
	data_file.close();
}
///////////////////////////////////////////write binary file
void writeBinaryForReco(char *OutputFileDir,int ProjAngle,char *InputNameFile,float v_strips[4][15000000],float t_strips[4][15000000],int energy[15000000][5],int GoodEventCounter)
{
char OutputFilename[512];
sprintf(OutputFilename, "%s/projection_%03.0f.bin",OutputFileDir,(double) ProjAngle);
////////////////////////////////////////////////////TV correction
	float par[5];
	CalF* f2cal[5];
	char TVfile[256];
	char CalibFile[256];
	sprintf(TVfile,"%s/TVcalib.txt",OutputFileDir);
	cout<<"Opening TV correction file: "<<TVfile<<endl;
	fstream TVcalfile(TVfile, std::ios_base::in); //text file with TV-corr data
	if(!TVcalfile.is_open() ){
		perror ("Error opening TV correction file");
		exit(1);
	}

	for (int i=0; i<5; ++i) {            // read fit parameters:
		TVcalfile >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] ;
	    f2cal[i]=new CalF(i,par);  // initialize t-v calibration  function
	}
	TVcalfile.close();
/////////////////////////////////////////// Prepare stuff for WEPL calibration
	float wepl_par[45];
	float Wet;
	//open text file with WEPL calibration data (parameters of 9 pol4 curves)
	sprintf(CalibFile,"%s/wet5calibExp.txt",OutputFileDir);
	fstream WEPLcalfile(CalibFile, std::ios_base::in); //sept
	cout<<"Opening WEPL calibration file: "<<CalibFile<<endl;
	if(!WEPLcalfile.is_open() ){
		perror ("Error opening WEPL calibration file");
		exit(1);
	}

  	for (int i=0; i<45; ++i) WEPLcalfile >> wepl_par[i];
 	WEPLcalfile.close();
	Wepl * WEPL=new Wepl(wepl_par); //initialize WEPL calibr.
	WEPL->SetEthresholds(3.,3.,3.,3.,3.); // Set all stage thresholds to 1 MeV

//////////////////////opening analysis file


//Histograms Initializtion
double Z=216.9 + 40; // approx position for the calorimeter entrance
//double ped[5] = {121.3, -71.5, -1137, 346.2, -49.}; // predefined pedestal values JULY2014 LLUMC
// double ped[5] = {431,-130,-20,224,60}; //SEPTEMBER 2014 LLUMC
double ped[5] = {549,92,204,575,385};// May and July 2015 CDH
float V[4], T[4]; 									// strip numbers or mm coordinates for the 4 t and V planes
float *V0, *T0, *u0;
float *V1, *T1, *u1;
float *V2, *T2, *u2;
float *V3, *T3, *u3;
float *WetBinary;
float u[4]={-214.35,-164.3,164.3,214.35};			// average coordinate of the T V boards for each plane
float Ped[5], Ene[5]; 								// pedestal and adc value to be converted in MeV and in WEPL
double Tcorr,Vcorr;									// Extrapolated coordinates used for TV correction
float Plane3[3], Plane4[3];							// Arrays used to plot the proton track between plane3 and plane4
													// and to intepolate the trajectory with a str8 line
int EvtNum=0;
int Event_Counter=GoodEventCounter;								//used to write data to binary files
V0= (float*) calloc (Event_Counter,sizeof(float));
V1= (float*) calloc (Event_Counter,sizeof(float));
V2= (float*) calloc (Event_Counter,sizeof(float));
V3= (float*) calloc (Event_Counter,sizeof(float));

T0= (float*) calloc (Event_Counter,sizeof(float));
T1= (float*) calloc (Event_Counter,sizeof(float));
T2= (float*) calloc (Event_Counter,sizeof(float));
T3= (float*) calloc (Event_Counter,sizeof(float));

u0= (float*) calloc (Event_Counter,sizeof(float));
u1= (float*) calloc (Event_Counter,sizeof(float));
u2= (float*) calloc (Event_Counter,sizeof(float));
u3= (float*) calloc (Event_Counter,sizeof(float));
WetBinary= (float*) calloc (Event_Counter,sizeof(float));

Plane3[2]=u[2];
Plane4[2]=u[3];

for (int EvtNum=0;EvtNum<GoodEventCounter;EvtNum++){
		if(EvtNum%100000==0)cout<<EvtNum<<endl;
		V0[EvtNum]= v_strips[0][EvtNum];
		V1[EvtNum]= v_strips[1][EvtNum];
		V2[EvtNum]= v_strips[2][EvtNum];
		V3[EvtNum]= v_strips[3][EvtNum];


		T0[EvtNum]=t_strips[0][EvtNum];
		T1[EvtNum]=t_strips[1][EvtNum];
		T2[EvtNum]=t_strips[2][EvtNum];
		T3[EvtNum]=t_strips[3][EvtNum];

		u0[EvtNum]=u[0];
		u1[EvtNum]=u[1];
		u2[EvtNum]=u[2];
		u3[EvtNum]=u[3];

		Plane3[0]=t_strips[2][EvtNum];
		Plane3[1]=v_strips[2][EvtNum];
		Plane4[0]=t_strips[3][EvtNum];
		Plane4[1]=v_strips[3][EvtNum];
		GetTVcoord(Z,Tcorr,Vcorr,Plane3,Plane4);//calling the function to project the beam

		for(int i=0;i<5;i++){
			Ene[i]=(float)energy[EvtNum][i];
			Ene[i]-=ped[i];//subtracting ped values
		   	Ene[i]=f2cal[i]->CalTVf(Tcorr,Vcorr,Ene[i]);//TVcorrection and MeV conversion
		}
		Wet=WEPL->EtoWEPL(Ene);
		WetBinary[EvtNum]=Wet;

	}
WriteBinaryFile ((float)ProjAngle, InputNameFile, Event_Counter, V0,V1,V2,V3, T0,T1,T2,T3, u0,u1,u2,u3, WetBinary,OutputFilename);

}
