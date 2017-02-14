%Function Process_SMOSxL1C to process SMOS L1C dual or four Stokes polarization Brightness Temperatures with multiple incidence angles. 
% USAGE 1:
% > [TSF, SSI] = Process_SMOSxL1C;
% then, a file browser will be open and select a DBL file to open.
% USAGE 2:
% > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS,outdir);
% WHERE:
% 'fname.DBL': full path string of DBL file to work with,
% GEOLIMS: a 4 element vector as [LAT1, LON1, LAT2, LON2],
% with elements indicating the latitude and longitude of bottom left
% of box and latitude and longitude of upper rigth of box.
% outdir: (optional) string with directory where filename.mat is stored. 
% USAGE 3:
% > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS);
% same as usage #2, but an output .mat file will not be generated.
% OUTPUTS:
% TSF: Structure with all variables for [Temp_Swath_Full] dataset.
% SSI: Structure with all variables for [Snapshot_Information]
% dataset (this output variable is optional).
% * Additionally, TSF includes the Brightness Temperatue in HV-pol after
% performing the Faraday rotation from the XY polarization plane.
% * The TSF memeber variables TB_Fixed_IncAngle and Fixed_IncAngle,
% are reprecenting the Brightness Temperatures at homogeneously spaced
% fixed incident angles (normally from 0 to 60 deg with 1 deg step).
% 
% More information: https://github.com/pablosaa/Process_SMOS_L1C
% Contact: Pablo Saavedra G., email: pablosaa@uni-bonn.de

