//*******************************************************************************************
//* New Program used to read raw data form binary file and store them in a useful format    *
//* P. Piersimoni August 05th, 2014 - Loma Linda University                                 *
//* Missing due to Si gaps 4th  hit coordinate reconstruction by 3-point quaratic spline    *
//* added by V.Bashkirov on June 28, 2015                                                   *
//* Latest review: P. Piersimoni Semptember 30th, 2015, UCSF, San Francisco, CA             *
//*******************************************************************************************

#include <iostream>
#include <stdlib.h>     // exit, EXIT_FAILURE
#include <fstream>
#include <string>
#include <bitset>       // std::bitset
#include <math.h>
#include <iomanip>      // std::setprecision
#include <stdio.h>
#include <cstdio>
#include <cerrno>
#include <vector>
#include <algorithm>
#include <limits.h>

#include <iterator>
#include <streambuf>
#include <sstream>
#include <functional>   // std::bind

// #include "TNtupleD.h"  // cern ROOT stuff
// #include <TCanvas.h>
// #include <TH2.h>
// #include <TFile.h>

#include "writeBinary.cc"
//*************************  if ROOT _NOT_ in compile using: g++ -std=c++0x -o codePP  pierReviewed_int.cpp
//************************** if ROOT in compile using: g++ -std=c++0x -c `root-config --cflags` PierVlad3points.cpp ; g++ `root-config --glibs` -o codePP PierVlad3points.o


using namespace std;
typedef unsigned char byte;
std::size_t BITS_PER_BYTE = CHAR_BIT;
//********************************************/global declaration
int stream_position=0;
unsigned int* buffer;
FILE* in_file;
bool PlaneFlag[4]={};
bool OneFPGAflag[12]={};
float v_strips[4][15000000];
float t_strips[4][15000000];
bool stripsDistance_flag=0;
int FPGAs_flag=0;
float StripDistance=5.;//
int checkWTHisgoingon=0;

int GoodEventCounter=0;
int Rec3EventCounter=0;
bool stop_reading=0;
size_t file_size=0;
int event_counter=0;
int pedestal[15000000][5];
int energy[15000000][5];
////////////////////////////////used for the 4th point reconstruction
float BoardT0[2]={-211.4,9999.};
float BoardT1[2]={-161.4,9999.};
float BoardT2[2]={161.4,9999.};
float BoardT3[2]={211.4,9999.};

float BoardV0[2]={-217.3,9999.};
float BoardV1[2]={-167.2,9999.};
float BoardV2[2]={167.2,9999.};
float BoardV3[2]={217.3,9999.};
//VB: define here distances between planes to speed up missing hit calculations.
double  z01, z02, z12, z13, z20, z21, z23;

////////////////////////////////////////////////////////////conversion to mm
float u[4]={-214.35,-164.3,164.3,214.35};//average coordinate of the T V boards for each plane
float tPlaneFirstStrip[4][4];//planes *sensors
float FirstStrip1[4]={176.644, 88.375, 0.089,  -88.182}; //first strip position for each sensor plane 1
float FirstStrip2[4]={172.77,  84.49,  -3.787, -92.074}; //first strip position for each sensor plane 2
float FirstStrip3[4]={187.643, 99.304, 11.028, -77.241}; //first strip position for each sensor plane 3
float FirstStrip4[4]={183.455, 95.195, 6.936,  -81.334}; //first strip position for each sensor plane 4
////////////////////////////////////////////////////////////conversion to mm
//********************************************/

////////union is used to reverse byte order
union int_reversal
{
    unsigned int int_val;
    struct
    {
        unsigned short int_hi;
        unsigned short int_lo;
    };
    char each_byte[4];//unused
};


unsigned short reverse_short_bytes(unsigned short x)
{
    //This function takes 16-bit short and switches the order of the high byte (first 8-bits) and low byte (last 8-bits)
	//shift low byte to high byte position
	//shift high byte to low byte position
	//recombine high/low bytes into 2 byte variable
    //return reversed_bytes;
    return ( x << 8 ) + ( x >> 8 );  //OK
}

unsigned int reverse_int_bytes(unsigned int x)
{
	//union
    int_reversal input_int;
    input_int.int_val = x;							// I put the x inside unsigned int (aa bb cc dd)
    unsigned short temp_short = input_int.int_hi;	// I take the first one half and I store it in a temp_short = aa bb
    input_int.int_hi = input_int.int_lo;			// I take the second one half and I put it in the first half int (cc dd cc dd)
    input_int.int_lo = temp_short;					// I put the temp_short in the second half	 int(cc dd aa bb)
    input_int.int_hi = reverse_short_bytes(input_int.int_hi); // I swap the first one half... (dd cc)
    input_int.int_lo = reverse_short_bytes(input_int.int_lo); //  ...and the second (bb aa)
    return input_int.int_val;            					  // I return the correct integer (dd cc bb aa)
}
                                             //current              //qued
void read_append_data( unsigned long long & bit_container, unsigned int & num_bits , int &origin)
{                             //this is 8 bits
    unsigned int shift_by = BITS_PER_BYTE * sizeof(*buffer); // shift is equivalent to the number of bits in a            byte (8) multiplied by the number of bytes
															 // composing the buffer (4)= 32 bits in total

    bit_container <<= shift_by;   					// Shift current_bits over by 4-bytes so another 4-bytes can be appended to the end

//	unsigned int*
	buffer = (unsigned int*) malloc (sizeof(unsigned int));

    fread (buffer,sizeof(*buffer),1,in_file);  		// Read next 4-bytes from the current reading position inside the file
//	fseek ( in_file, sizeof(*buffer), origin );
	origin= origin+sizeof(*buffer);


    auto temp = reverse_int_bytes(buffer[0]);       // Reverse the order of the bytes previously read into the buffer
    bit_container |= temp;                      	// Add these bits to the end of current_bits into the 4-bytes slot
													// opened up by shifting the existing bits by 4-bytes
    num_bits += shift_by;                        	// Having just read and appended 4-bytes = 32-bits to current_bits,
													// add 32 to the count of # bits represented by current_bits
	if(origin>file_size )stop_reading=1;
	free(buffer);
}
                                             //current              //qued
unsigned int extract_N_bits( unsigned long long & bit_container, unsigned int & num_bits, unsigned int bits_2_extract  )
{
    unsigned long long extracted_bits = bit_container;
	if (num_bits<bits_2_extract)read_append_data( bit_container, num_bits , stream_position);
    unsigned int shift_size = num_bits - bits_2_extract;
    extracted_bits >>= shift_size;
    bit_container = bit_container - ( extracted_bits << shift_size );
    num_bits -= bits_2_extract;
    return extracted_bits;
}
//this function groups together the read and extract functions, and makes the check on the queud bits
unsigned int just_read(unsigned long long & bit_container, unsigned int & num_bits, unsigned int bits_2_extract, int &origin)
{
	if(origin>file_size ){
		stop_reading=1;
		return 0;
		}else
		if(num_bits<bits_2_extract)read_append_data( bit_container, num_bits, origin );
		unsigned int required_bts = extract_N_bits(bit_container, num_bits, bits_2_extract);
		return required_bts;

}


//////////////////////////////////////////////////////////file header function
void read_file_header(unsigned long long fileHeader_bits, char *namefile)
{
//the number 13784398 corresponds to the bit string: 0000 0000 1101 0010 0101 0101 0100 1110 = "1RUN" in std ASCII( in a fancy way)
	if(fileHeader_bits== 13784398)cout<<"Reading pCT file: "<<namefile<<endl;
	if(fileHeader_bits!= 13784398){
		perror ("The file is not a pCT file");
		exit(1);
	}
}
//////////////////////////////////////////////////////////run number function
void read_run_number(unsigned long long runNumber_bits){cout<<"Run Number: "<<runNumber_bits<<endl;}
//////////////////////////////////////////////////////////start time function
void read_run_startTime(unsigned long long startTime_bits)
{
//old version to be used for sept and july 2014 data
//	union {
//			float f;
//			unsigned int integer;
//		  } float2int;
//	float2int.integer = startTime_bits;
//	int TIME = float2int.f;
	time_t RealTime =startTime_bits;
	cout<<"The run start time is: "<< ctime (&RealTime)<<endl;
}

//////////////////////////////////////////////////////////status bits function
void read_statusBits(unsigned long long status_bits){if (status_bits==1)cout<<"Time tag = 1"<<endl;}

//////////////////////////////////////////////////////////program version function
void read_programVersion(unsigned long long programVersion_bits){cout<<"The program version is: "<<programVersion_bits<<endl;}
//////////////////////////////////////////////////////////projection angle
void read_projection_angle(unsigned long long projectionangle_bits){cout<<"The projection angle is: "<<projectionangle_bits <<endl;}



//////////////////////////////////////////////////////////file header function
// I need a lot of parameters because I have to look for the string.
// I keep in memory the old values of current and queued and if I don't find the STRING("1pCT") I restore the old values,
// but I sottract 1 to queued so that at the next step (I prevent a new reading) bits will be extract starting after 1 position:
// if I have in memory e.g. 0000>0000STRING ( > is the point marked by queued), I will extract the first time 0000STRING.
// if I sottract 1 to queued at the next step, after restoring current, I will have: 00000>000STRING, so I will extract 000STRING.
// I keep doing this way until I extract exactly STRING&mask, because in the container probaly there no anly zerosw before queued_bits
// NB I need to restore current and queued because the extract function automatically delete the extract bits  from current
// and update queued

 																				 //current                                   //queued
bool read_BegOfEvent (unsigned long long BegOfEvent_bits, unsigned long long & bit_container, unsigned long long temp_cont , unsigned int & num_bits, unsigned int temp_queu)
{
//the string is: 1111 0000 0100 0011 0101 0100= "1pCT" in ASCII std= 15745876
	bool Found=0;
	unsigned long long mask =pow(2,24)-1;
	if((BegOfEvent_bits & mask) == 15745876){
		event_counter++;
 		Found=1;
	}else{
		bit_container= temp_cont ;
		num_bits=temp_queu-4;//I can do 4-bit steps, because the BoE usually starts at the byte or half byte beginning
	}
	return Found;
}

//////////////////////////////////////////////////////////time tag function
void read_timeTag(unsigned long long timeTag_bits) {//cout<<"Event time tag (in clock cycles): "<< timeTag_bits  <<endl
;	}

//////////////////////////////////////////////////////////time since previous trigger function
void read_timeDelta(unsigned long long timeDelta_bits) {//cout<<"Time since the previous trigger (in clock cycles): "<< timeDelta_bits  <<endl
;	}

//////////////////////////////////////////////////////////event header function
// the event header is made of 6 pieces (pz):
//	-start bits : 10 (11+22zeros)
//	-4 error flags (1+21,20,19,18, zeros)
//	-18 bits: event counter.(18 ones)
// I use a mask for each piece to extract the desired bits and I shift  right the extracted bits to have an integer

void read_eventHeader (unsigned long long eventHeader_bits)
{
	unsigned long long eventHeader_mask[6]={0xC00000,0x200000,0x100000,0x80000,0x40000,0x3FFFF};// rapresentaion HEX to do a mask
	unsigned long long eventHeader_pz[6]={};
	for (int l=0;l<6; l++)eventHeader_pz[l]=eventHeader_bits & eventHeader_mask[l];

	if((eventHeader_pz[1] >>= 21 )==1 )cout<<"***********Error flag: !Incorrect FPGA adress received!***"<<endl;
	if((eventHeader_pz[2] >>= 20 )==1 )cout<<"***********Error flag: !Tag mismatch error!***"<<endl;
	if((eventHeader_pz[3] >>= 19 )==1 )cout<<"***********Error flag: !CRC error!***"<<endl;
	if((eventHeader_pz[4] >>= 18 )==1 )cout<<"***********Error flag: !Chip error!***"<<endl;
	if((eventHeader_pz[0] >>= 22 )==2&&(eventHeader_pz[5]%100000)==0 ) cout<<"Starting event id "<<eventHeader_pz[5]<<" counter "<< event_counter<< endl;}
//////////////////////////////////////////////////////////FPGAs header function (12 bits)
// FPGAs header is made of 6 pieces (pz):
//	-4 bits : FPGA adress 0-11, 1111+8zeros
//	-4 error flags (1+7,6,5,4, zeros){ignored}
//	-4 bits: cips number.(1111)
// I use a mask for each piece to extract the desired bits and I shift  right the extracted bits to have an integer
int read_FPGAsHeader (unsigned long long FPGAsHeader_bits, int FPGA_num)
{
	unsigned long long FPGAsHeader_mask[6]={0xF00,0x80,0x40,0x20,0x10,0xF};	// rapresentaion HEX to do a mask
																		   	// here I am ignoring some pieces
	unsigned long long FPGAsHeader_pz[6]={};
	for (int l=0;l<6; l++)FPGAsHeader_pz[l]=FPGAsHeader_bits & FPGAsHeader_mask[l];
	if (FPGA_num == (FPGAsHeader_pz[0] >>= 8)){
//		cout<<"reading FPGA "<<FPGAsHeader_pz[0]<<"..."<<endl;
		return FPGAsHeader_pz[5];
	}else{
		perror ("Error: FPGA Address mismatch");
		exit(1);
	}

}

// ASICs header is made of 6 pieces (pz):
//	-2 bits (overflow+unused):  11+10zeros=0xC00
//	-4 bits: number of clusters (111+6zeros)=0x3C0
//	-2bits:Chip error and parity (ignored)=0x30
//	-4 bits: chips number.(1111)=0xF
// I use a mask for each piece to extract the desired bits and I shift  right the extracted bits to have an integer
int read_ASICHeader (unsigned long long ASICHeader_bits, int &ASIC_add)
{
	unsigned long long ASICHeader_mask[4]={0xC00,0x3C0,0x30,0xF};// representaion HEX to do a mask
	unsigned long long ASICHeader_pz[4]={};
	for (int l=0;l<4; l++)ASICHeader_pz[l]=ASICHeader_bits & ASICHeader_mask[l];
	ASIC_add= ASICHeader_pz[3] ;
	return (ASICHeader_pz[1] >>= 6);
}

// strip header is made of 2 pieces (pz):
//	-6 bits :  111111+6zeros = 0xFC0
//	-6 bits :  111111 = 0x3F
// I use a mask for each piece to extract the desired bits and I shift  right the extracted bits to have an integer
int read_stripHeader (unsigned long long stripHeader_bits, int &strip_add)
{
	unsigned long long stripHeader_mask[2]={0xFC0,0x3F};// representaion HEX to do a mask
	unsigned long long stripHeader_pz[2]={};
	for (int l=0;l<2; l++)stripHeader_pz[l]=stripHeader_bits & stripHeader_mask[l];
	strip_add= stripHeader_pz[1] ;
	return (stripHeader_pz[0] >>= 6);
}
// Energy FPGAs header is made of 6 pieces (pz):
//	-1 bit : pedestal_flag( 1= ped, 0=no ped) 1+12zeros=0x800
//	-4 bits :  front-ent tag 1111+7zeros=0x780
// 	-1 bit: data type ( 0= samples, 1= reduced) 1+6zeros= 0x40
//	-1 bit: error flag 1+5zeros =0x20
//	-2 bits: number of channel per fpga (3 or 2) 11+000= 0x18
//	-3 bits: OTR( not used if reduced is active)111=7
// I use a mask for each piece to extract the desired bits and I shift  right the extracted bits to have an integer
int read_EnFPGAs (unsigned long long EnFPGAs_bits, bool &ped_flag)
{
	unsigned long long EnFPGAs_mask[6]={0x800,0x780,0x40,0x20,0x18,7};// representaion HEX to do a mask
	unsigned long long EnFPGAs_pz[6]={};
	for (int l=0;l<6; l++)EnFPGAs_pz[l]=EnFPGAs_bits & EnFPGAs_mask[l];
	ped_flag = (EnFPGAs_pz[0] >>= 11);
	int frontend_tag= (EnFPGAs_pz[1]>>= 7);
	int Entype_flag =(EnFPGAs_pz[2]>>= 6);
	bool EnFPGAs_err=(EnFPGAs_pz[3]>>= 5);
	if(EnFPGAs_err==1)	perror ("Error flag in Energy FPGA ");
	int channels_num=(EnFPGAs_pz[4]>>= 3);
	int OTR=(EnFPGAs_pz[5]);

	return	channels_num;

}
float ConvertTomm(float StripT,int plane)// ussed for the convertion to mm only for strip T, the conversion for V is trivial
{
		if (StripT >=0   && StripT<384 )StripT=tPlaneFirstStrip[plane][0]-.228*StripT;
		if (StripT>=384  && StripT<768 )StripT=tPlaneFirstStrip[plane][1]-.228*(StripT-384);
		if (StripT>=768  && StripT<1152)StripT=tPlaneFirstStrip[plane][2]-.228*(StripT-768);
		if (StripT>=1152 && StripT<1536)StripT=tPlaneFirstStrip[plane][3]-.228*(StripT-1152);

	return StripT;

}


// this is a function calculating the correct strip number to be kept: among all the hit strips, I determine Max and min and
// I calculate the distance between max and min, if this distance is <S tripDistance set above, I set the Stripflag to 1. If two chips are hit, may be
// they are contigous and so may be just 1 proton hit them, but if the distance is too high 2 protons hit them.
// I can say this just from the strip, I don't need a control on the chip numbers.
// Obviousely if the hit strip is just one this fuction return it, but I don't need another if inside the code
int keep_correct_strip(vector<int> &strip_numbVec,size_t size_vec)
{

 	stripsDistance_flag=0;
	float Max, min;
	Max=min=strip_numbVec[0];
	int true_strip;
	for(int i=0; i<size_vec; i++) {
  		if( strip_numbVec[i]>Max ) Max=strip_numbVec[i];
  		if( strip_numbVec[i]<min ) min=strip_numbVec[i];
	}
	if ((Max-min)<StripDistance) stripsDistance_flag=1;
	true_strip=lround((Max+min)/2);
	return true_strip;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is used to fill the correct strip array with the calculated strip number, according to FPGA adress after

void fill_strips_arrays(float strip_calc, int FPGANo)
{
    if(FPGANo<4){//V boards
        // converting V strip no to mm
		strip_calc=-43.717+strip_calc*.228;
		v_strips[FPGANo][GoodEventCounter]=strip_calc;
		if(FPGANo==0)BoardV0[1]=strip_calc;
		if(FPGANo==1)BoardV1[1]=strip_calc;
		if(FPGANo==2)BoardV2[1]=strip_calc;
		if(FPGANo==3)BoardV3[1]=strip_calc;
	}
    //T boards
	else if(FPGANo==4  || FPGANo==5) {
        strip_calc=ConvertTomm(strip_calc,0);
		t_strips[0][GoodEventCounter]=strip_calc;
		BoardT0[1]=strip_calc;
	}else if(FPGANo==6  || FPGANo==7) {
		strip_calc=ConvertTomm(strip_calc,1);
		t_strips[1][GoodEventCounter]=strip_calc;
		BoardT1[1]=strip_calc;
	}else if(FPGANo==8  || FPGANo==9)  {
		strip_calc=ConvertTomm(strip_calc,2);
		t_strips[2][GoodEventCounter]=strip_calc;
		BoardT2[1]=strip_calc;
	}else if(FPGANo==10 || FPGANo==11) {
		strip_calc=ConvertTomm(strip_calc,3);
		t_strips[3][GoodEventCounter]=strip_calc;
		BoardT3[1]=strip_calc;
    }
}
void Set_FPGAsFlag(int &FPGAs_flag, int FPGANo)//used to check how many hits I have
{

	OneFPGAflag[FPGANo]=1;
	PlaneFlag[0]= OneFPGAflag[0] && (OneFPGAflag[4]!=OneFPGAflag[5]);
	PlaneFlag[1]= OneFPGAflag[1] && (OneFPGAflag[6]!=OneFPGAflag[7]);
	PlaneFlag[2]= OneFPGAflag[2] && (OneFPGAflag[8]!=OneFPGAflag[9]);
	PlaneFlag[3]= OneFPGAflag[3] && (OneFPGAflag[10]!=OneFPGAflag[11]);
	FPGAs_flag= PlaneFlag[0]+PlaneFlag[1]+PlaneFlag[2]+PlaneFlag[3];
}

void reset_FPGAflags()
{
    for (int i=0;i<12;i++)OneFPGAflag[i]=0;
    for (int i=0;i<4;i++)PlaneFlag[i]=0;
    BoardT0[1]=9999.;
    BoardT1[1]=9999.;
    BoardT2[1]=9999.;
    BoardT3[1]=9999.;
    BoardV0[1]=9999.;
    BoardV1[1]=9999.;
    BoardV2[1]=9999.;
    BoardV3[1]=9999.;
}


int GetForthPoint()
{
	double x0, x1, x2, x3;
	double z0, z1, z2, z3;
	double a,b,c;
	if(BoardT0[1]==9999.){
        z0=BoardT0[0];
    	z1=BoardT1[0];
    	z2=BoardT2[0];
    	x1=BoardT1[1];
    	x2=BoardT2[1];
    	x3=BoardT3[1];
    	a = ((x2 - x3) *z23-(x2 - x1)*z21) *z21;
    	b = (x2 - x3) *z23 - 2 * a*z2;
    	BoardT0[1] = x1+(2*a*z1+b)/z01;
     	t_strips[0][GoodEventCounter]=BoardT0[1]; ;
    // for histogramm
    	return 0;
	} else if(BoardV0[1]==9999.){
  	    z1=BoardV1[0];
		z2=BoardV2[0];
 		x1=BoardV1[1];
		x2=BoardV2[1];
		x3=BoardV3[1];
		a = ((x2 - x3) *z23-(x2 - x1)*z21) *z21;
		b = (x2 - x3) *z23 - 2 * a*z2;
		BoardV0[1]=x1+(2*a*z1+b)/z01;
		v_strips[0][GoodEventCounter]=BoardV0[1]; ;
	    return 10;
	} else if(BoardT1[1]==9999.){
        z1=BoardT1[0];
        z2=BoardT2[0];
        x0=BoardT0[1];
        x2=BoardT2[1];
        x3=BoardT3[1];
        a = ((x2 - x3)*z23-(x2 - x0)*z20 ) *z20;
        b = (x2 - x3)*z23 - 2 * a*z2;
        c = x2 - a*z2*z2 - b*z2;
        BoardT1[1] = a*z1*z1 + b*z1 + c;
        t_strips[1][GoodEventCounter]=BoardT1[1]; ;
	       return 1;
	} else if(BoardV1[1]==9999.){
 		z1= BoardV1[0];
		z2=BoardV2[0];
		x0=BoardV0[1];
		x2=BoardV2[1];
		x3=BoardV3[1];
		a = ((x2 - x3)*z23-(x2 - x0)*z20 ) *z20;
		b = (x2 - x3)*z23 - 2 * a*z2;
		c = x2 - a*z2*z2 - b*z2;
		BoardV1[1]=a*z1*z1 + b*z1 + c;
 		v_strips[1][GoodEventCounter]=BoardV1[1];
	return 11;
    } else if(BoardT2[1]==9999.){
 		z1=BoardT1[0];
		z2=BoardT2[0];
		x0=BoardT0[1];
		x1=BoardT1[1];
		x3=BoardT3[1];
		a = ((x0 - x1) *z01-(x1 - x3)*z13) *z13;
		b = (x0 - x1) *z01 - 2 * a*z1;
		c = x1 - a*z1*z1 - b*z1;
		BoardT2[1]=a*z2*z2 + b*z2 + c;
		t_strips[2][GoodEventCounter]=BoardT2[1];
	return 2;
    } else if(BoardV2[1]==9999.){
 		z1=BoardV1[0];
		z2= BoardV2[0];
 		x0=BoardV0[1];
		x1=BoardV1[1];
		x3=BoardV3[1];
		a = ((x0 - x1) *z01-(x1 - x3)*z13) *z13;
		b = (x0 - x1) *z01 - 2 * a*z1;
		c = x1 - a*z1*z1 - b*z1;
		BoardV2[1]=a*z2*z2 + b*z2 + c;
 		v_strips[2][GoodEventCounter]=BoardV2[1];
        // FPGAs_flag=4;
		return 12;
    } else if(BoardT3[1]==9999.){
 		z1=BoardT1[0];
		z2=BoardT2[0];
		x0=BoardT0[1];
		x1=BoardT1[1];
		x2=BoardT2[1];
		a = ((x0 - x1)*z01-(x1 - x2)*z12) *z12;
		b = (x0 - x1) *z01 - 2 * a*z1;
		BoardT3[1]=x2-(2*a*z2+b)/z23;
 		t_strips[3][GoodEventCounter]=BoardT3[1];
        // FPGAs_flag=4;
		return 3;
	} else if(BoardV3[1]==9999.){
		z1=BoardV1[0];
		z2=BoardV2[0];
		x0=BoardV0[1];
		x1=BoardV1[1];
		x2=BoardV2[1];
		a = ((x0 - x1)*z01-(x1 - x2)*z12) *z12;
		b = (x0 - x1) *z01 - 2 * a*z1;
		BoardV3[1]=x2-(2*a*z2+b)/z23;
 		v_strips[3][GoodEventCounter]=BoardV3[1];
        // FPGAs_flag=4;
		return 13;
	}  else return 9999;
}
bool checkCosAnddist()
{
    bool cutsFlag =0;
    float P0[3]={t_strips[0][GoodEventCounter],v_strips[0][GoodEventCounter],u[0]};
    float P1[3]={t_strips[1][GoodEventCounter],v_strips[1][GoodEventCounter],u[1]};
    float P2[3]={t_strips[2][GoodEventCounter],v_strips[2][GoodEventCounter],u[2]};
    float P3[3]={t_strips[3][GoodEventCounter],v_strips[3][GoodEventCounter],u[3]};
    //check if the cos(angle) betwwen the two straight lines is greater than 0.98 and if the distance between the straight line at 0 is less than 10mm
    double vectorU[3];
    double vectorV[3];
    for (int i=0;i<3;i++){
      vectorU[i]=P1[i]-P0[i];
      vectorV[i]=P3[i]-P2[i];
    }
    double UdotV=vectorU[0]*vectorV[0]+vectorU[1]*vectorV[1]+vectorU[2]*vectorV[2];
    double modU=sqrt(vectorU[0]*vectorU[0]+vectorU[1]*vectorU[1]+vectorU[2]*vectorU[2]);
    double modV=sqrt(vectorV[0]*vectorV[0]+vectorV[1]*vectorV[1]+vectorV[2]*vectorV[2]);
    double cosTheta = UdotV/(modU*modV);
    double PinAt0[2]={};
    double PoutAt0[2]={};
    GetTVcoord(0.,PinAt0[0],PinAt0[1],P1,P0);
    GetTVcoord(0.,PoutAt0[0],PoutAt0[1],P3,P2);
    double distanceAt0= sqrt( (PinAt0[0]-PoutAt0[0])*(PinAt0[0]-PoutAt0[0]) +(PinAt0[1]-PoutAt0[1])*(PinAt0[1]-PoutAt0[1]) );
    if(cosTheta>=0.98 && distanceAt0<10. ){
        cutsFlag =1;
    }
    return cutsFlag;
}

int Get_StripNumber(int ChipNo, int ChannelNo, int FPGANo ){
//The strip number depends on the chip Number(0-11)and the Channel Number(0-63).
//Each FPGA has 12 (0-11) chips(asics) each of them contains 64(0-63) channels
//T planes are read by 2 FPGAs, V planes by 1.
//Because we want each plane to give us integer numbers (strip#) with this formula we got:
//For V planes: values 0-383 top to bottom
//For T planes: values 0-1535 left to right
//NOTE that for T planes the orientation changes from front and rear planes (mirrored)

//OUTPUT by Andriy
// V: 0..383  FPGA#3    384..767
// T: 767---#11--0 767---#10---0
// V: 0..383  FPGA#2   384..767
// T: 767---#9--0  767---#8----0
//
// T: 0---#6---767 0----#7----767
// V: 767..384  FPGA#1    383..0
// T: 0----#4--767 0----#5----767
// V: 767..384  FPGA#0    383..0

//What we got
// ^
// | V3: 0..383  FPGA#3    0..383
// | T: 0---#11--767 767---#10---1535    		*****REAR FPGAs
// | V2: 0..383  FPGA#2   0....383
// | T: 0---#9--767  768---#8----1535
// |
// | T: 0---#6---767 768----#7----1535    same
// | V1: 0....383  FPGA#1    0...383			*****Front FPGAs
// | T: 0----#4--767 768---#5--1535   L to R
// | V0: 0...383  FPGA#0    0...383   Top to bot
//beam direction
	int StripNo=0;
	if(FPGANo>=0&&FPGANo<=3){
	//Vplanes
    	if(ChipNo<=5&&ChipNo>=0)StripNo= (5-ChipNo)*64+ChannelNo; 		//Left	 |for all the 4 V planes
    	if(ChipNo<=11&&ChipNo>=6)StripNo=(ChipNo-6)*64+(64-ChannelNo-1);//Right	 |is the same
	}
	if(FPGANo>=4&&FPGANo<=7){
//	T Planes Front, before the the detector
    	if(FPGANo%2==0)	StripNo= ChipNo*64+(64-ChannelNo-1);//left
    	if(FPGANo%2!=0)	StripNo= 768+ ChipNo*64+(64-ChannelNo-1);//right
	}
	if(FPGANo>=8&&FPGANo<=11){
//	T Planes rear, after the the detector
    	if(FPGANo%2==0)	StripNo=768+ (11-ChipNo)*64+ChannelNo;//right
    	if(FPGANo%2!=0)	StripNo=(11-ChipNo)*64+ChannelNo ;    //left
	}
	return StripNo;

}

void write_data2file(string dirOut, string NameFile)
{
	dirOut = dirOut+"/TextFiles";
	string command = "mkdir -p "+dirOut;
	system(command.c_str());
	ofstream output;
	NameFile= dirOut+"/"+ NameFile+".out";
	output.open(NameFile.c_str(),ios::trunc);
	if (output.is_open()){
		cout <<"Writing data to file: "<<NameFile<<endl;
	}else {
		cout <<"Can't write data to file: "<<NameFile<<endl;
		perror ("Error opening file");
		exit(1);
	}
	output<<"EvtNb\t V0mm \t T0mm \t V1mm \t T1mm \t V2mm \t T2mm \t V3mm \t T3mm \t Ped0 \t ADC0 \t Ped1 \t ADC1 \t Ped2 \t ADC2 \t Ped3 \t ADC3  \t Ped4 \t  ADC4"<<endl;
    output<<GoodEventCounter<<endl;
	for(int i=0;i<GoodEventCounter;i++){
		output<<i;
		for (int j=0;j<4;j++)output<<"\t"<< v_strips[j][i]<< "\t"<<t_strips[j][i];
		for (int j=0;j<5;j++)output<<"\t"<< pedestal[i][j]<<"\t"<<energy[i][j];
		output<<endl;
	}
	output.close();
}

//*          |\\  //|    //\    ||  |\\  ||  * * * *__________ *************************************************
//*          ||\\//||   // \\   ||  ||\\ ||   * * * __________ *************************************************
//*          || \/ ||  //ZZZ\\  ||  || \\||  * * * *__________ *************************************************
//*          ||    || //     \\ ||  ||  \\|  _________________ *************************************************
//*****************************************  _________________ *************************************************

int main (int argc, char* argv[])
{
///////////////////////////////////Needed for conversion to mm
for (int i=0; i<4; ++i){
tPlaneFirstStrip[0][i]=FirstStrip1[i];//just one array for all the planes
tPlaneFirstStrip[1][i]=FirstStrip2[i];
tPlaneFirstStrip[2][i]=FirstStrip3[i];
tPlaneFirstStrip[3][i]=FirstStrip4[i];
}
///////////////////////////////////Needed to speed-up missing hit reconstruction
z01 = 1 / (BoardT0[0] - BoardT1[0]);
z02 = 1 / (BoardT0[0] - BoardT2[0]);
z12 = 1 / (BoardT1[0] - BoardT2[0]);
z13 = 1 / (BoardT1[0] - BoardT3[0]);
z23 = 1 /(BoardT2[0] - BoardT3[0]);
z20 = -z02;
z21 = -z12;
int Pnum;
///////////////////**************************
///// Root file and Histograms booking-Used to test the 4th point recon:DON'T Delete it
// TFile *fcn = new TFile("Recon3p.root","RECREATE");
// 	cout<< "Opened TFile Recon3p.root" <<endl;
// TH1F* hT[4];
// TH1F* hV[4];
// for (int i=0; i<4; ++i) {
//     hT[i] = new TH1F(Form("hT%d",i), Form("T reconstructed for plane %d",i), 800,-200,200 );
//     hV[i] = new TH1F(Form("hV%d",i), Form("V reconstructed for plane %d",i), 200,-50,50 );
// 	}
///////////////////**************************

////opening the binary file
bool binary=0;
char data_filename[256];
char completeNameFile[256];
char Inputdir[500];
char Outputdir[500];

for (int i=0; i<argc; i++)printf("%s\n",argv[i]);
if (argc<6) {
cout<< "Please add parameters to the executable:\n " << argv[0] << " projection# fileName[Object_RUN#]  RawDataDir[path/to/MMDDYYYY/RUN#] date[MMDDYYYY] OutputMode[0 =ASCII 1=BINARY]" <<endl;
return 0;
}

int AngleNb = atoi(argv[1]);
binary = atoi(argv[5]);
cout << "Binary boolean is " << binary << endl;
cout<<"Projection n "<<AngleNb<<endl;
sprintf( data_filename, "%s_%03.0f",argv[2],(double)AngleNb);
sprintf(Inputdir,"%s/Input",argv[3]);
sprintf( completeNameFile, "%s/%s.dat",Inputdir,data_filename);
sprintf(Outputdir,"%s/Output/%s",argv[3],argv[4]);
cout<<"Reading "<<completeNameFile <<endl;
in_file = fopen(completeNameFile, "rb");
if (in_file == NULL){
	perror ("Error opening file");
	exit(1);
}
fseek(in_file, 0L, SEEK_END);
file_size = ftell(in_file);
rewind(in_file);
//cout<<"size"<<file_size<<endl;
unsigned long long current_bits = 0; 	// holds the ULL representation of all the bits read into memory
										// and placed in the correct order
unsigned int queued_bits = 0;			// # of bits currently represented by the ULL variable "current_bits"
										//(i.e. # bits in queue waiting to be processed)
unsigned long long extracted_bits;      // bit containing information needed
int required_bits=0;
//this for extracts the first 8 byte from file
for( int j = 0; j < 2; j++){
read_append_data( current_bits, queued_bits, stream_position );
//	cout<<j<<"j "<<"current "<<(bitset<64>)current_bits<<endl;
}

//Reading FILE header, calling read_file_header function
//file header is made by 4 bytes identifien a pCT file, 32 bits needed
required_bits=32;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_file_header(extracted_bits,data_filename);

//reading Run Number: 24 bits
required_bits=24;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_run_number(extracted_bits);

//reading Run start time:32 bits
required_bits=32;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_run_startTime(extracted_bits);

//reading status bits:8 bits
required_bits=8;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_statusBits(extracted_bits);

//reading program version number:8 bits
required_bits=8;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_programVersion(extracted_bits);
//reading projection_angle:12 bits
required_bits=12;
extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
read_projection_angle(extracted_bits);

//////////////////////////////////////////////////////////////////////////////////////////here a cycle all over the while loop all over the events
while(!stop_reading){

//reading Beginning-of-Event identifier: 24 bits. the string is: 1111 0000 0100 0011 0101 0100= "1pCT" in ASCII std= 15745876(int dec)
//I need a for cycle because in the documentation it is not said how many 0 bits are inserted before the "1pCT" string
//
	bool Eureka=0;
	required_bits=24;
	while(!Eureka){
	if(stream_position>=file_size )break;
		if(queued_bits<required_bits)read_append_data( current_bits, queued_bits, stream_position );
		unsigned long long temp_container= current_bits;
		unsigned int temp_queued = queued_bits;
		extracted_bits = extract_N_bits( current_bits, queued_bits, required_bits );
		Eureka = read_BegOfEvent(extracted_bits,current_bits, temp_container,queued_bits, temp_queued);
		}
	if(stream_position>=file_size )break;
//reading the  event time tag: 36 bits
	required_bits=36;
	extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
	read_timeTag(extracted_bits);

//reading the time since previous trigger: 12 bits
	required_bits=12;
	extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
	read_timeDelta (extracted_bits);
//reading Event header: 24 bits
	required_bits=24;
	extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
	read_eventHeader (extracted_bits);

//reading 12 tracker FPGA headers: 12 bits

////////////////////////////
// these variables are used to determine wheter the event is a "good one" or not.
// FPGAS flag says if the proton hits 4 planes and excludes 2 fpga in the same plane
// strip_numb_vec is used to store multiple strip hit
// stripsDistance_flag tells if the distance between two hit strips (determined by StripDistance) is too much for a single proton
// if both FPGAs and strip are true the event is declred good
//		bool FPGAs_flag=0;
	FPGAs_flag=0;
	vector<int> strip_numb_vec;
//////////////////////////////////
	int Strip=0;
	int hit_chips=0;
	int strip_address=0;
	int clusters_num=0;
	int strip_num=0;
	int FPGA_numb;
	int ASIC_address=0;
	required_bits=12;
	for (FPGA_numb=0;FPGA_numb<12;FPGA_numb++){
		    size_t strip_counter=0;
			strip_numb_vec.clear();
        extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
			hit_chips=	read_FPGAsHeader (extracted_bits, FPGA_numb);
			if (hit_chips>0){
//this cycle reads chip-headers, according to FPGANofChips: 12 bits.
				required_bits=12;
				for(int i= 0;i<hit_chips;i++){
        		extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
				clusters_num = read_ASICHeader (extracted_bits, ASIC_address);
				required_bits=12;
                //this cycle reads Strips-headers, according to FPGANofStrips.
				while(clusters_num>0){//
						extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
						strip_num =  read_stripHeader (extracted_bits, strip_address);//if strip number is 0 means 1 strip is hit
					for(int i=0;i<=strip_num;i++){
							Strip=Get_StripNumber(ASIC_address, strip_address+i,FPGA_numb);
							strip_numb_vec.push_back(Strip);
							strip_counter++;
					}
  					clusters_num--;
  				}
  			}
  		}
  		//function to determine which strip to be kept
  		strip_numb_vec.resize(strip_counter);
  		if(strip_counter!=0){
  			Strip=keep_correct_strip(strip_numb_vec, strip_counter);
            if(stripsDistance_flag==1)fill_strips_arrays((float)Strip, FPGA_numb);
  		} else stripsDistance_flag=0;
        if(stripsDistance_flag==1)Set_FPGAsFlag(FPGAs_flag, FPGA_numb);
    }//ending the FPGA reading
//////////////////////////////////////////////////////
//2 energy FPGAs reading
	int count_channel=0;
	for(int i=0;i<2;i++){
	bool pedestal_flag=0;
	required_bits=12;
		extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
	int channels=read_EnFPGAs(extracted_bits, pedestal_flag);
		for (int j=0;j<channels;j++){
		if (pedestal_flag){
			required_bits=8;
			extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
			int value =(int8_t)extracted_bits;
//				cout<<"  ped "<<(bitset<8>)extracted_bits<<" v "<<value<<endl;
			pedestal[GoodEventCounter][count_channel]=value;
    //			cout<<"  ped array "<<pedestal[GoodEventCounter][count_channel]<<endl;
		}
		required_bits=16;
		extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
		int value1 =(int16_t)extracted_bits;
//			cout<<"  ene "<<(bitset<16>)extracted_bits<<" v "<<value1<< endl;
		energy[GoodEventCounter][count_channel]=value1;
  //			cout<<"  ene array "<<energy[GoodEventCounter][count_channel]<<endl;
		count_channel++;
		if (pedestal_flag==0 && i==1){
			required_bits=4;
			extracted_bits =just_read(current_bits, queued_bits,required_bits, stream_position);
		}
	}
}
if(FPGAs_flag==3) {
  Rec3EventCounter++;
  Pnum=GetForthPoint();
  if(FPGAs_flag==3) {
    Rec3EventCounter++;
    Pnum=GetForthPoint();
    if (Pnum!=9999){
       //  	if(Pnum>9) hV[Pnum-10]->Fill(v_strips[Pnum-10][GoodEventCounter]);
       //  	else hT[Pnum]->Fill(t_strips[Pnum][GoodEventCounter]);
          if(checkCosAnddist()){
              FPGAs_flag=4;
              checkWTHisgoingon++;
          }
      }
  }

if(FPGAs_flag==4)GoodEventCounter++;
reset_FPGAflags();
//////////////////////////////////////////////////////
//	if(stream_position%10000==0)cout<<stream_position<<" "<<file_size<<endl<<endl;
if(stream_position>file_size )stop_reading=1;
if(GoodEventCounter%100000==0)cout<<GoodEventCounter<<" GoodEventCounter"<<" Stream position "<<stream_position<<" file size  "<<file_size<<endl<<endl;
}

fclose(in_file);
if(!binary)write_data2file(Outputdir,data_filename);
if(binary){
	writeBinaryForReco(Outputdir, AngleNb, completeNameFile,v_strips,t_strips,energy,GoodEventCounter);
}
// for (int ical=0; ical<4; ++ical) {
// new TCanvas;
// hT[ical]->Draw(" ");  hV[ical]->Draw(" ");}
// fcn->Write(); fcn->Close("R");
cout<<event_counter<<" events processed; "<<GoodEventCounter<<" good events stored;" <<endl;
cout<<"Efficency "<<((float)GoodEventCounter)/(float)event_counter*100.<<"\%"<< endl;
cout<<"3points found "<<(float)(Rec3EventCounter)/(float)event_counter*100.<<"\%"<< endl;
cout<<"4 point found but discarded  "<<(float)(Rec3EventCounter- checkWTHisgoingon)/(float)event_counter*100.<<"\%"<< endl;
cout<<"Exit program, everything went good, Good Bye! ;)"<<endl;
return 0;
}
