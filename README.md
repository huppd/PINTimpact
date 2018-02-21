PImpact is a Periodic Incompressible navier--stokes solver on Massively PArallel
CompuTers.

folder structure:

- run: collection of running scripts
- src: scr code
  - src_c: c++ code 
  - src_f: Fortran code
  - test: collection of many unit tests
- XML: colection of different parameters and cases

 files types
 - f90: fortran code
 - hpp: c++ header files
 - cpp: c++ source files
 - cxx: c++ main files
 - xml: parameter xml files
 - py: python scripts
 - txt: text files for cmake


 installation:

 - needed tools: cmake + compilers
 - needed libraries: blas/lapack, mpi, hdf5
 - trilinos needed packages(Teuchos, Belos, NOX)

