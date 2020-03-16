#!/bin/csh
NULIB_LOC = /Users/jspark971/external/NuLib
HDF5_LOC = /Users/jspark971/external/hdf5-1.10.5/hdf5
# export DYLD_LIBRARY_PATH=/Users/jspark971/Documents/research/isotropicSQA/IsotropicSQA/external/hdf5-1.10.5/hdf5/lib/

INCLUDE  = -I${HDF5_LOC}/include

LIBRARY = ${NULIB_LOC}/src/nulib.a -L${HDF5_LOC}/lib -lhdf5 -lhdf5_fortran -lhdf5_cpp -L/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9 -lgfortran -lomp

WHICHCODE = sqa2.psma

COMP = g++ -g -std=gnu++11 -O3 -Wall -Xpreprocessor -fopenmp # -DNDEBUG

${WHICHCODE}.o: ${WHICHCODE}.cpp
	rm -f sqa2.psma.x sqa.tar.gz
	tar --exclude=*.lum --exclude=*.cyl -cvzf sqa.tar.gz ./*
	${COMP} ${WHICHCODE}.cpp -o ${WHICHCODE}.x ${INCLUDE} ${LIBRARY}
