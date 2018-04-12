**Pimpact** stands for  <b>p</b>eriodic <b>i</b>ncompressible Navier--Stokes solvers on
<b>m</b>assively <b>pa</b>llel <b>c</b>ompu<b>t</b>ers.


## 1) Obtaining code

To get the code, one has to execute

<code> git clone https://github.com/huppd/PINTimpact.git </code>

in a shell.
This will download the folder *PINTimpact* into the folder of execution.
This folder has the following folder structure

- *run*<!--: collection of running scripts-->
- *src*<!--: scr code-->
  - *src_c*<!--: c++ code -->
  - *src_f*<!--: Fortran code-->
  - *test*<!--: collection of many unit tests-->
- *XML*<!--: collection of different parameters and cases-->

The *run* folder contains various python scripts, that can start multiple jobs and create
folder hierarchies for parameter studies or scaling tests.  The *src* folder contains
three subfolders, *src_c* for the C++ part of the code, *src\_f* for the Fortran90 part of
the code, and *test* for the unit-tests.  The *XML* folder contains various *xml*
parameter files.  These are used to set up the parameters and problems for the *Pimpact*
solver.

<!--### Files types-->
 <!--- *f90*: fortran code-->
 <!--- *hpp*: c++ header files-->
 <!--- *cpp*: c++ source files-->
 <!--- *cxx*: c++ main files-->
 <!--- *xml*: parameter xml files-->
 <!--- *py*: python scripts-->
 <!--- *txt*: text files for cmake-->

### Generating documentation

The documentation can be created by executing

<code> doxygen Doxyfile </code>

in the *PINTimpact* folder.
This will generate the documentation in the *doc* folder.
The documentation can be read by opening the file *doc/html/index.html* by a
preferred browser.
Further information about compiling and using the code can be found there.



## 2) Compiling code

We recommend to compile the code in two versions: a release version that is optimized, and
a debug version that provided runtime tests and further debugging informations.  A user,
that is only running the code without adjustments might, do without a debug version.  A
developer, who is improving the performance of the code, might even want to use a third
profiling version.

*Pimpact* is written partly in C++ and Fortran90, so we recommend to use a compiler that
can compile both and directly link those.  Test have been done with the intel compiler and
the gcc compiler.

The following libraries are needed:
 - needed tools: cmake + compilers
 - needed libraries: blas/lapack, mpi, hdf5
 - trilinos needed packages(Teuchos, Belos, NOX)

### 2a) Installing *Trilinos*

Depending if your platform already provides 

### 2b) Compiling *Pimpact*

<!--debug/release-->

## 3) Running tests 

## 4) Running code

## 5) Running parameter studies

## 6) Analyzing and Visualizing
