# Process_SMOS_L1C
Octave/Matlab rapit processor for SMOS satellite brightness temperature L1C data

Copyright (c) 2016-2017 Pablo Saavedra G.  
version 1.0, see LICENSE  

## COMPILATION:  ##
For compilation the GNU Scientific Library is needed. In case of Octave the MATIO library is additionally needed.
### from linux console:  ###

      >[MATLAB_BIN_PATH]/mex CFLAGS='$CFLAGS -std=gnu++11' Process_SMOSxL1C.cpp -lgsl -lgslcblas  
creates the mex function `Process_SMOSxL1C.mexa64`  

### from MATLAB workspace:  ###

      > mex Process_SMOSxL1C.cpp -lgsl -lgslcblas
creates the mex function `Process_SMOSxL1C.mexa64`  

### from OCTAVE workspace:  ###

      > mkfileobj Process_SMOSxL1C.cpp -lgsl -lgslcblas -lmatio  
this creates the octave mex `function Process_SMOSxL1c.mex`  
### Makefile option ###
Use the provided `Makefile` adapting the enviroment variable MATLABROOT or OCTAVEROOT with the path for matlab version installed in your system, e.g. MATLABROOT = /usr/local/matlab7.8/, then in linux console run:

      > make
to create the matlab mex function `Process_SMOSxL1C.mexa64` or

      > make octave
to create the Octave mex function `Process_SMOSxL1C.mex`

## HOW TO USE ##
### USAGE 1:###

      > [TSF, SSI] = Process_SMOSxL1C;  
then, a file browser will be open and select a DBL file to open.  
### USAGE 2:  ###

      > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS,outdir);  
WHERE:  
'fname.DBL': full path string of DBL file to work with,  
GEOLIMS: a 4 element vector as [LAT1, LON1, LAT2, LON2],  
with elements indicating the latitude and longitude of bottom left of box and latitude and longitude of upper rigth of box.  
outdir: (optional) string with directory where filename.mat is stored.  
### USAGE 3:  ###

      > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS);  
same as usage #2, but an output .mat file will not be generated.  
OUTPUTS:  
TSF: Structure with all variables for [Temp_Swath_Full] dataset.  
SSI: Structure with all variables for [Snapshot_Information] dataset (this output variable is optional).  
* Additionally, TSF includes the Brightness Temperatue in HV-pol after performing the Faraday rotation from the XY polarization plane.  
* The TSF memeber variables TB_Fixed_IncAngle and Fixed_IncAngle, are reprecenting the Brightness Temperatures at homogeneously spaced fixed incident angles (normally from 0 to 65 deg with 1 deg step).  

## EXAMPLES ##

      > TSF=Process_SMOSxL1C('/mission/smos/SM_OPER_MIR_SCLF1C_20150702T042618_20150702T051937_620_001_1.DBL',[10 -5 60 20],'/data/output/');  
The example above process the DataBlock file SM_OPER_MIR_SCLF1C_20150702T042618_20150702T051937_620_001_1.DBL located at /mission/smos/ directory, extracting data from a geographically limited box with 10° latitude and -5° longitude left-bottom boundary and 60° latitude and 20° longitude rigth-upper boundary. Then the processed data is archived at the /data/ouput directory as a MAT-file, and delivered to workspace in the TSF structure variable with fields given as foolow:    

      TSF = 

                  File_Name: 'SM_OPER_MIR_SCLF1C_20150702T042618_20150702T051937_620_001_1.DBL'
               GEOBoxLimits: [2x2 double]
             Start_End_Time: [2x17 char]
               GridPoint_ID: [1041x1 double]
         GridPoint_Latitude: [1041x1 double]
        GridPoint_Longitude: [1041x1 double]
         GridPoint_Altitude: [1041x1 double]
             GridPoint_Mask: [1041x1 double]
            BT_Data_Counter: [1041x1 double]
                      flags: {1041x1 cell}
                 BTvalue_re: {1041x1 cell}
                 BTvalue_im: {1041x1 cell}
        PixelRadiometry_acc: {1041x1 cell}
            Incidence_angle: {1041x1 cell}
              Azimuth_angle: {1041x1 cell}
      FaradayRotation_angle: {1041x1 cell}
    GeometricRotation_angle: {1041x1 cell}
         Snapshot_ID_Pixels: {1041x1 cell}
            Footprint_Axis1: {1041x1 cell}
            Footprint_Axis2: {1041x1 cell}
                       TBxy: {1041x1 cell}
                       TBhv: {1041x1 cell}
                         RA: {1041x1 cell}
             idx_SnapshotID: {1041x1 cell}
          TB_Fixed_IncAngle: [1041x61x4 double]
             Fixed_IncAngle: [61x1 double]

