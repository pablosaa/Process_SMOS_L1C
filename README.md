# Process_SMOS_L1C
<<<<<<< HEAD
Octave/Matlab rapid processor for SMOS satellite brightness temperature L1C data
=======
Octave/Matlab rapit processor for SMOS satellite full- or dual-polirized brightness temperature L1C data
>>>>>>> 1b71c170faf1ab3eaa4fba0a91c66fc0f68e892f


Process_SMOSxL1C is a MEX-function written in C++ to process L1C scientific data producs from the MIRAS interferometer on board of the Soil Moisture Ocean Salinity (SMOS) satellite, a European Space Agency mission.
The L1C data product comes with data for brightness temperature at L-band in dual-polarization or full-polarization (4 Stokes vector) mode. 

The MEX-function comes in two flavors: as MATLAB function and as a GNU/Octave function. The MATLAB function has as file name "Process\_SMOSxL1C.mexa64", while the Octave function is named as "Process\_SMOSxL1C.mex".

The Function has been tested under MATLAB version 2007, 2010b, 2016a. As well as Octave versions 3.6.2 and 4.0.3.

For information about how to compile and examples how to run the function and plot content, see the [wiki](https://github.com/pablosaa/Process_SMOS_L1C/wiki) page in this repository.

This software is registered at [Zenodo.org](https://zenodo.org) under [![DOI](https://zenodo.org/badge/81903236.svg)](https://zenodo.org/badge/latestdoi/81903236)

For technical details about what is implemented in the MEX-function, see the Open Access paper:

**Saavedra, P. & Simmer, C.**, (2018). _An Octave/MATLAB® Interface for Rapid Processing of SMOS L1C Full Polarization Brightness Temperature_. Journal of Open Research Software. 6(1), p.2. [DOI:10.5334/jors.165](http://doi.org/10.5334/jors.165).

Copyright (c) 2016-2017 Pablo Saavedra G.  
version 1.0, see LICENSE  
