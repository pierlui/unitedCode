# unitedCode
code for pCT raw data reading
Please add parameters to the executable:
 ./a.out (executable)
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
