# Makefile to compile the MEX function <Process_SMOSxL1C>  from linux console:
#// Copyright (c) 2016 Pablo Saavedra Garfias
#//
#// this file is part of 'Process_SMOSxL1C', see LICENSE.

HOST = $(shell hostname)

ifeq ($(HOST),mysu.yolle)
	MATLABROOT = /usr/local/MATLAB/R2010b
else
	MATLABROOT = /usr/local/matlab7.8
#	MATLABROOT = /usr/local/matlab9.1
endif

# General Options to pass to the compiler defined in CFLAGS
COMMON_CXX =  -O2 -Wall -Wextra -std=gnu++11
MYLIBS = -lgsl -lgslcblas

# Matlab Options:
MYFLAGS = CXXFLAGS='$$CXXFLAGS $$COMMON_CXX' #-Wall -O2 -std=gnu++11'
CC = $(MATLABROOT)/bin/mex # Matlab compiler

# Octave Options:
OCTAVEROOT = /usr/bin
MYFLAGS_OCT = -DOCT_MEX_FILE=1 --mex
CC_OCT = $(OCTAVEROOT)/mkoctfile  # Octave compiler
OCT_CXX = `mkoctfile -p CXXFLAGS` $(COMMON_CXX)
MYOCTLIBS = $(MYLIBS) -lmatio
SOURCE = Process_SMOSxL1C.cpp

# Rules to build MATLAB MEX function:
all: Process_SMOSxL1C.mexa64

Process_SMOSxL1C.mexa64: $(SOURCE) 
	 $(CC) $(MYFLAGS) $^ $(MYLIBS)

# Rules to build Octave MEX function:
octave: Process_SMOSxL1C.mex

Process_SMOSxL1C.mex: $(SOURCE)
	CXXFLAGS="$(OCT_CXX)" \
	$(CC_OCT) $(MYFLAGS_OCT) $^ $(MYOCTLIBS)


clean:
	rm -f Process_SMOSxL1C.o Process_SMOSxL1C.mex Process_SMOSxL1C.mexa64
