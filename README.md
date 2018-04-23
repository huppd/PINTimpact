**Pimpact** stands for  <b>p</b>eriodic <b>i</b>ncompressible Navier--Stokes solvers on
<b>m</b>assively <b>pa</b>llel <b>c</b>ompu<b>t</b>ers.


## 1) Obtaining code

To get the code, one has to execute

``` git clone https://github.com/huppd/PINTimpact.git ```

in a terminal.  This will download the folder *PINTimpact* into the folder of execution.
This folder has the following folder structure

- *run*
- *src*
- *src_c*
- *src_f*
- *test*
- *XML*

The *run* folder contains various python scripts, that can start multiple jobs and create
folder hierarchies for parameter studies or scaling tests.  The *src* folder contains
three subfolders, *src\_c* for the C++ part of the code, *src\_f* for the Fortran90 part
of the code, and *test* for the unit-tests.  The *XML* folder contains various *xml*
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

Depending if your platform already provides an actuall verision of  *Trilinos* (12.1)
including the packages *Teuchos*, *Belos*, *NOX*, this step can be skipped and directly
proceeded with **2b)**.  A quick installation guide can be found at
https://github.com/trilinos/Trilinos/blob/master/INSTALL.rst.  A full installation
reference can be found at https://trilinos.org/docs/files/TrilinosBuildReference.html.  We
provide here short customized guide for the needs of Pimpact.
  
The Trilinos code can be downloaded by 

``` git clone https://github.com/trilinos/Trilinos.git ```

We reccomend to save the following steps in a bash script for each platform.  To install
Trilinos after downloading the code, three steps have to be done: make the buildsytem by
`cmake`, build the code by `make`, install the build by `make install`.  Optionally the
build can be tested by `ctest`.
  
``` cmake $ARGS ${TRILINOS_SRC_HOME} ```

The `TRILINOS_SRC_HOME` is the path to the downloaded folder. If the git clone comand from
above has been executed in the homefolder, then it is `TRILINOS_SRC_HOME="~/Trilinos"`.


  ``` -D CMAKE_INSTALL_PREFIX="~/trilinos/debug"  ```

  ``` -D CMAKE_BUILD_TYPE:STRING=DEBUG \ ```

  ``` -D Trilinos_ENABLE_DEBUG=ON \ ```

  ``` -D TPL_ENABLE_MPI:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_TESTS:BOOL=ON \ ```

  ``` -D BUILD_SHARED_LIBS:BOOL=OFF \ ```

  ``` -D TPL_ENABLE_LAPACK:BOOL=ON \ ```

  ``` -D TPL_ENABLE_Boost:BOOL=OFF \ ```

  ``` -D TPL_ENABLE_BLAS:BOOL=ON \ ```

  ``` -D TPL_ENABLE_HDF5:BOOL=OFF \ ```

  ``` -D TPL_ENABLE_MKL:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_CHECKED_STL:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_CXX11:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_OpenMP:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Anasazi:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_AztecOO:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Fortran:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_Teuchos:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_Kokkos:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Belos:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_NOX:BOOL=ON \ ```

  ``` -D Trilinos_ENABLE_Epetra:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_PyTrilinos=OFF \ ```

  ``` -D Trilinos_ENABLE_ML:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Thyra:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Tpetra:BOOL=OFF \ ```

  ``` -D Trilinos_ENABLE_Galeri:BOOL=OFF \ ```


  After a successful creating the build system, Trilinos can be compiled.

  ``` make -j4 ```

  Optional: if the `Trilinos_ENABLE_TESTS:BOOL=ON` has been set the Trilinos build can be
  tested in the temporary build director by

  ``` ctest . ```

  After the succesfull build, Trilinos can be installed by

  ``` make install ```

  This will basically move the necessary header and library files to the set
  `CMAKE_INSTALL_PREFIX` folder.  The temporary Trilions folder can be removed or cleaned
  afterwards.


### 2b) Compiling *Pimpact*

  We discurage in builts, as cmake is generating plenty of temporal files and folders.
  We recomend having one or multiple build director in the *PINTimpact* folder, called
  accordingly debug, release, or profile.

  This comad has to be executed in the build folder for example in `~/PINTimpact/release`.
  The last argument for `cmake` has to be the path to the main `CMakeLists.txt`, which can
  be for example `~/PINTimpact/src`

  ``` cmake src_path ```

  ``` -DTrilinos_DIR=... ```

  ``` -DfindHDF5=ON/OFF ```

  ``` -DHDF5_INCLUDE_DIRS=... ```

  ``` -DHDF5_LIBRARY_DIRS=... ```

  ``` -DCMAKE_BUILD_TYPE=Release ```

  ``` -DCMAKE_BUILD_TYPE=Debug ```

  ``` -DCMAKE_BUILD_TYPE=Profile ```


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


## 5) Running parameter studies

    In the folder *run* are multiple python scripts to run parameter studies.



## 6) Analyzing and Visualizing

  Depending on the set parameters in the xml file different output files are produced.
  There are different txt files that contain informations about the performance of the
  solver such as achived tolerance or number of iterations.
  Different `h5`files can be wirtten out, for restarting, analyzing, or visualizing. 
  If the fields are not written for restarting they are interpolated on the pressure grid
  and are accompanied with `xmf` files.
  These `xmf` files can be opened with for example [paraview](https://www.paraview.org/)
  or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).

