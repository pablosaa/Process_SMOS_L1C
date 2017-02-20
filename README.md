# Process_SMOS_L1C
Octave/Matlab rapit processor for SMOS satellite brightness temperature L1C data


Process_SMOSxL1C is a MEX function written in C++ to process L1C scientific data producs from the MIRAS interferometer on board of the Soil Moisture Ocean Salinity (SMOS) satellite, a European Space Agency mission.
The L1C data product comes with data for brightness temperature at L-band in dual-polarization or full-polarization (4 Stokes vector) mode. 

The MEX-function comes in two flavors: as MATLAB function and as a GNU Octave function. The MATLAB function has file name of "Process_SMOSxL1C.mexa64", while the Octave function is named as "Process_SMOSxL1C.mex".

The Function has been tested under MATLAB version 2007, 2010b, 2016a.

For information about how to compile and examples how to run the function, see the wiki page in this repository.

Copyright (c) 2016-2017 Pablo Saavedra G.  
version 1.0, see LICENSE  
