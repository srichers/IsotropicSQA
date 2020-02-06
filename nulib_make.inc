F90=mpif90
F90FLAGS= -O0 -fopenmp -g

#if you want to take advantage of openmp, set OPENMP=1
OPENMP=1

#if you want to take advantage of mpi, set MPI=1
#mpif90 must be in your PATH (openmpi)
MPI = 0

MODINC="-I ./"

#You must have a HDF5 version installed with
#the _same_ compiler, this often means compiling
#it by yourself.
#See README file for compiling HDF5

#HDF5DIR=/usr/local/hdf5-1.8.17_gnu-5.4.0/hdf5
HDF5INCS=-I$(HDF5_DIR)/include
HDF5LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz -Wl,-rpath=$(HDF5_DIR)/lib

#additional flag needed to include other's work, like Hempel's EOS or Timmes' Helmholtz EOS
NUCLEI_HEMPEL = 0
HELMHOLTZ_EOS = 0
WEAK_RATES = 0
DEBUG = 0