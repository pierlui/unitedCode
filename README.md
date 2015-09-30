# unitedCode
code for pCT raw data reading
 
The directory where the program is rum must contain the 2 classes for WEPL reconstruction and TV (Wepl.h, CalTV.h) and the writeBinary.cc file.
compile the main code by: g++ -std=c++0x  -o codePP maincode.cpp.
At the moment two main codes are uploaded in this repository:
1: pierReview_int.cpp
   this program just keeps the proton history in which 4 good hits are registered
2: PreProcessing3point.cc
   this program recovers the 4th hit coordinates when 3 good hits are registered

Please add parameters to the executable:
 .codePP (executable)
 projection# 
 fileName[Object_RUN#]  
 RawDataDir[path/to/MMDDYYYY/RUN#] 
 date[MMDDYYYY] 
 OutputMode[0 =ASCII 1=BINARY]
 
This program is meant to be used in the pCT file storage system, i.e.:
the input file dir is supposed to be:
path/to/object/DATE/RUN#/Input
the input file name is:
Object_Run#_projectionAngle.dat
The Output directory is:
path/to/object/DATE/RUN#/Output/DATE/
In case the OuputMode is set to ASCII a new Dir containing the output text file will be created e.g.,:
path/to/object/DATE/RUN#/Output/DATE/TextFiles/Object_Run#_projectionAngle.out

The Output/DATE directory must contain the WEPL calibration and TVcorrection coefficient files:
path/to/object/DATE/RUN#/Output/DATE/wet5calibExp.txt
path/to/object/DATE/RUN#/Output/DATE/TVcalib.txt

The pedestal values 

