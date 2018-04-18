**Pimpact** stands for  <b>p</b>eriodic <b>i</b>ncompressible Navier--Stokes solvers on
<b>m</b>assively <b>pa</b>llel <b>c</b>ompu<b>t</b>ers.


## 1) Obtaining code

To get the code, one has to execute

``` git clone https://github.com/huppd/PINTimpact.git ```

in a shell.  This will download the folder *PINTimpact* into the folder of execution.
This folder has the following folder structure

- *run*
- *src*
- *src_c*
- *src_f*
- *test*
- *XML*

The *run* folder contains various python scripts, that can start multiple jobs and create
folder hierarchies for parameter studies or scaling tests.  The *src* folder contains
three subfolders, *src_c* for the C++ part of the code, *src\_f* for the Fortran90 part of
the code, and *test* for the unit-tests.  The *XML* folder contains various *xml*
parameter files.  These are used to set up the parameters and problems for the *Pimpact*
solver.

### File types
- *f90*: fortran code
- *hpp*: c++ header files
- *cpp*: c++ source files
- *cxx*: c++ main files
- *xml*: parameter xml files
- *py*: python scripts
- *txt*: text files for cmake

### Generating documentation

The documentation can be created by executing

``` doxygen Doxyfile ``` 

in the *PINTimpact* folder.
This will generate the documentation in the *doc* folder.
The documentation can be read by opening the file *doc/html/index.html* by a
preferred browser.
Further information about compiling and using the code can be found there.
This has to be redone, if changes in the code documentation have been done.


## 2) Compiling code

We recommend to compile the code in two versions: a release version that is optimized, and
a debug version that provides runtime tests and further debugging informations.  A user,
  that is only running the code without adjustments might do without a debug version.  A
  developer, who is working on the performance of the code, might even want to use a third
  version for profiling.

  *Pimpact* is written partly in C++ and Fortran90, so we recommend to use a compiler that
  can compile both and directly link those.  Test have been done with the intel compiler and
  the gcc compiler.

  Needed libraries are MPI, HDF5, LAPACK, and trilinos.

  We are using functionality from the ... Standard so, we recomend using OpenMPI ... or
  MPICH ....

  The recomended HDF5 version is 1.8.13


### 2a) Installing *Trilinos*

  Depending if your platform already provides 

- trilinos needed packages(Teuchos, Belos, NOX)

### 2b) Compiling *Pimpact*

  We discurage in builts, as cmake is generating plenty of temporal files and folders.
  We recomend having one or multiple build director in the *PINTimpact* folder, called
  accordingly debug, release, or profile.


## 3) Running tests 

  After finishing step **2b)** the implementation can be tested. One has to change into the
  build director and can simply run

  ``` ctest . ```

  This will run all test with multiple parameters.  The source code for the test can be find
  in *PINTImpact/src/test* folder.  They are written using the [unit testing
  support](https://trilinos.org/docs/dev/packages/teuchos/doc/html/group__Teuchos__UnitTest__grp.html)
  provided by the *Teuchos* package from *Trilinos*.
  .
  The according executable can be found in the build director */test*.  If the executable
  are run with a `--help` flag, they print all possible possible parameter for the command
  line.


## 4) Running code

  In the folder *PINTimpact* you can find different executable.
  The according binarires have the same names just without the suffix .cxx.
  To print the their comand line parameters they can be executed with the `--help` flag.
  This is an overview over the executables 

  - analyzer3D.cxx analyzes results from peri_navier3D
  - `--filename`
  - `--restart`
  - convDiff.cxx solver for convection-diffusion problems (has no comandline parameters)
  - convDiff2.cxx solver for convection-diffusion problems (has no comandline parameters)
  - convDiffMG.cxx solver for convection-diffusion problems (has no comandline parameters)
- convDiffMG2.cxx solver for convection-diffusion problems (has no comandline parameters)
  - modeConvDiff.cxx solver for harmonic convection-diffusion problems
  - `--filename`
  - `--realCase`
  - peri_navier2D.cxx spectral in time solver for two-dimensional Navier--Stokes equations
  - peri_navier3D.cxx spectral in time solver for the three-dimesnional Navier--Stokes equations
  - `--filename`
  - `--restart`
  - `--nf_restart`
  - peri_navier3DTime.cxx solver with finite differences in time
  - `--filename`
  - peri_navier4.cxx deprecated version of peri_navier3DTime.cxx
  - `--re`
  - `--alpha2`
  - `--rad`
  - `--amp`
  - `--xm`
  - `--ym`
  - `--flow`
  - `--forcing`
  - `--domain`
  - `--l1`
  - `--l2`
  - `--l3`
  - `--n1`
  - `--n2`
  - `--n3`
  - `--nt`
  - `--np1`
  - `--np2`
  - `--np3`
  - `--npt`
  - `--linSolName`
  - `--nonLinSolName`
  - `--lineSearchName`
  - `--fixType`
  - `--leftPrec`
  - `--precType`
  - `--maxIter`
  - `--tolBelos`
  - `--tolNox`
  - `--tolSchur`
  - `--tolPrec`


## 5) Running parameter studies


## 6) Analyzing and Visualizing

  Depending on the set parameters in the xml file different output files are produced.
  There are different txt files that contain informations about the 
