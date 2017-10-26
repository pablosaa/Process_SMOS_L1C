# Process_SMOS_L1C
Octave/Matlab rapit processor for SMOS satellite full- or dual-polirized brightness temperature L1C data


Process_SMOSxL1C is a MEX function written in C++ to process L1C scientific data producs from the MIRAS interferometer on board of the Soil Moisture Ocean Salinity (SMOS) satellite, a European Space Agency mission.
The L1C data product comes with data for brightness temperature at L-band in dual-polarization or full-polarization (4 Stokes vector) mode. 

The MEX-function comes in two flavors: as MATLAB function and as a GNU Octave function. The MATLAB function has file name of "Process_SMOSxL1C.mexa64", while the Octave function is named as "Process_SMOSxL1C.mex".

The Function has been tested under MATLAB version 2007, 2010b, 2016a. As well as Octave versions 3.6.2 and 4.0.3.

For information about how to compile and examples how to run the function and plot content, see the [wiki](https://github.com/pablosaa/Process_SMOS_L1C/wiki) page in this repository.

This software is registered at [Zenodo.org](https://zenodo.org) under [![DOI](https://zenodo.org/badge/81903236.svg)](https://zenodo.org/badge/latestdoi/81903236)

Copyright (c) 2016-2017 Pablo Saavedra G.  
version 1.0, see LICENSE  
