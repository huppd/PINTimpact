#*************************************************************************************************************
#* IMPACT                                                                                                    *
#* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
#* Mai 2005 - Dec 2011                                                                                       *
#*************************************************************************************************************


#--- compilation of the source code ---
#
# use this file to define the following macros for the compilation of the source code:
#    COMP (compiler)
#    DST  (target directory)
#    PRE  (preprocessor options)
#    CUSR (compile only usr-files ("FALSE" or anything else))
#    INCL (include options)
#    LIBS (library options)
#    OPT1 (compiler options for linking)
#    OPT2 (compiler options for compilation)
#
# note: - please edit only this file but not the file "compile"
#       - compilation and linking requires the external libraries MPI-2, HDF5 1.6 or 1.8 and Lapack
#       - normally you have to use the same compiler which was used for compilation of MPI-2, HDF5 and Lapack
#         on your computer (otherwise the compiler may not be able to correctly compile and link the code)
#       - the calls to the external libraries are written for double precision arithmetic (8 byte) by default
#         such that you have to modify them correspondingly (e.g. by means of a preprocessor macro, cf. the
#         example "xt5r4" below) if you want to compile the code for single precision
#       - for dynamic memory allocation you have to use the macro "PRE = -DALLOC"
#       - for non-Boussinesq simulations you have to use the macro "PRE = -DNONBOUSSINESQ" (or "PRE =
#         '-DNONBOUSSINESQ -DALLOC'" in case of dynamic memory allocation)
#       - to compile/link/clean only the usr-files use the macro "CUSR = FALSE"
#       - the source code uses free-length lines which may conflict with the default settings of your compiler
#         (e.g. gfortran complains), i.e. you have to set appropriate flags to avoid that
#       - to compile the code, type "gmake" or "gmake build"
#       - to delete the modules and objects in the target directory type "gmake clean"
#       - if you use more than one computer and/or compiler you can also work with targets (see the commented
#         out examples below and the man pages for "make" or "gmake", respectively)
#
DST = ../prog
CUSR = FALSE


#*************************************************************************************************************
hpc2:  COMP = mpif90
hpc2:  DST  = ../prog
hpc2:  #PRE  = -DALLOC
hpc2:  PRE  = -DALLOC
hpc2:  CUSR = FALSE
#hpc2:  INCL = -I/usr/include 
#hpc2:  INCL = -I/usr/include -I/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/gcc_4.4.6/include/ 
hpc2:  INCL = -I/usr/include  -I/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.6.2/gcc_4.7.2/include
#hpc2:  LIBS = -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
#hpc2:  LIBS = -L/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/gcc_4.4.6/lib -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
hpc2:  LIBS = -L/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.6.2/gcc_4.7.2/lib -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
hpc2:  OPT1 =  -ffree-line-length-none \
			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
			  -J$(DST)
hpc2:  OPT2 = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
			  -J$(DST) -xf95-cpp-input
hpc2:  build

hpc2usr:  COMP = mpif90
hpc2usr:  DST  = ../prog
hpc2usr:  #PRE  = -DALLOC
hpc2usr:  PRE  = -DALLOC
hpc2usr:  CUSR = TRUE
hpc2usr:  INCL = -I/usr/include 
hpc2usr:  LIBS = -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
hpc2usr:  OPT1 = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
			  -J$(DST)
hpc2usr:  OPT2 = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
			  -J$(DST) -xf95-cpp-input
hpc2usr:  build
#*************************************************************************************************************


#*************************************************************************************************************
brutus:  COMP = mpif90
#brutus:  COMP = h5pfc
brutus:  DST  = ../prog
brutus:  #PRE  = -DALLOC
brutus:  PRE  = -DALLOC
brutus:  CUSR = FALSE
#brutus:  INCL = -I/usr/include -I/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/pgi_12.6/include/ 
brutus:  INCL = -I/usr/include -I/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/gcc_4.4.6/include/ 
#brutus:  LIBS = -L/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/pgi_12.6/lib -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
brutus:  LIBS = -L/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/gcc_4.4.6/lib -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
brutus:  OPT1 = -r4 -i4 \
			  -O4 -fast -fastsse -Mipa=fast,inline  -Minline=reshape -Mvect=nosizelimit -Munroll -Mprefetch=n:80 -Mprefetch=w -Mflushz -Minform=inform \
			  -module $(DST) -Mpreprocess -Ktrap=divz,inv,ovf -Mdaz \
			  -DMPI_REAL8=MPI_REAL4 -DH5T_NATIVE_DOUBLE=H5T_NATIVE_REAL \
			  -DDGBTRF=SGBTRF -DDGBTRS=SGBTRS -DDSWAP=SSWAP -DDSCAL=SSCAL -DDGER=SGER -DDCOPY=SCOPY -DDLASWP=SLASWP -DDTRSM=STRSM -DDGEMM=SGEMM -DDCOPY=SCOPY
brutus:  OPT2 = $(OPT1)
brutus:  build
#*************************************************************************************************************


#*************************************************************************************************************
brutususr:  COMP = mpif90
#brutususr:  COMP = h5pfc
brutususr:  DST  = ../prog
brutususr:  #PRE  = -DALLOC
brutususr:  PRE  = -DALLOC
brutususr:  CUSR = TRUE
brutususr:  INCL = -I/usr/include -I/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/pgi_12.6/include/ 
brutususr:  LIBS = -L/cluster/apps/hdf5/1.8.8/x86_64/openmpi_1.4.5/pgi_12.6/lib -lhdf5_fortran -lhdf5 -lhdf5_hl  -lz -lm -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
brutususr:  OPT1 = -r4 -i4 \
			  -O4 -fast -fastsse -Mipa=fast,inline  -Minline=reshape -Mvect=nosizelimit -Munroll -Mprefetch=n:80 -Mprefetch=w -Mflushz -Minform=inform \
			  -module $(DST) -Mpreprocess -Ktrap=divz,inv,ovf -Mdaz \
			  -DMPI_REAL8=MPI_REAL4 -DH5T_NATIVE_DOUBLE=H5T_NATIVE_REAL \
			  -DDGBTRF=SGBTRF -DDGBTRS=SGBTRS -DDSWAP=SSWAP -DDSCAL=SSCAL -DDGER=SGER -DDCOPY=SCOPY -DDLASWP=SLASWP -DDTRSM=STRSM -DDGEMM=SGEMM -DDCOPY=SCOPY
brutususr:  OPT2 = $(OPT1)
brutususr:  build
#*************************************************************************************************************

##--- compile the code using targets (alternative) ---
## note: - this is an example for a Makefile with four different targets named "gfort", "xt5", "xt5r4", "p5"
##       - create your own targets or edit the examples
##       - define the macros, set the compilation target "build" and include the file "compile"
##       - to compile the code, type "gmake [your target]"
##-------------------------------------------------------------------------------------------------------------
#gfort:  COMP = gfortran
##gfort:  DST  = ../prog
#gfort:  DST = ./prog
#gfort:  PRE  = -DALLOC
#gfort:  CUSR = FALSE
#gfort:  INCL = -I/usr/include -I/usr/lib64/gfortran/modules
#gfort:  LIBS = -lmpi_f77 -lmpi -L/usr/lib64 -lhdf5_fortran -lhdf5 -lhdf5_hl -lz -lm -llapack
#gfort:  OPT1 = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
#			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
#			  -J$(DST)
#gfort:  OPT2 = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
#			  -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#			  -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno \
#			  -J$(DST) -xf95-cpp-input
#gfort:  build
##-------------------------------------------------------------------------------------------------------------
#xt5:    COMP = ftn
#xt5:    DST  = ../prog
#xt5:    PRE  = -DALLOC
#xt5:    CUSR = FALSE
#xt5:    INCL = -I/opt/cray/hdf5-parallel/1.8.3.0/hdf5-parallel-pgi/include
#xt5:    LIBS = -L/opt/cray/hdf5-parallel/1.8.3.0/hdf5-parallel-pgi/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm
#xt5:    OPT1 = -r8 -i4 \
#			  -O4 -fast -fastsse -Mipa=fast,inline -tp istanbul-64 -Minline=reshape -Mvect=nosizelimit -Munroll -Mprefetch=n:80 -Mprefetch=w -Mflushz -Minform=inform \
#			  -module $(DST) -Mpreprocess -Ktrap=divz,inv,ovf -Mdaz
#xt5:    OPT2 = $(OPT1)
#xt5:    build
##-------------------------------------------------------------------------------------------------------------
#xt5r4:  COMP = ftn
#xt5r4:  DST  = ../prog
#xt5r4:  PRE  = -DALLOC
#xt5r4:  CUSR = FALSE
#xt5r4:  INCL = -I/opt/cray/hdf5-parallel/1.8.3.0/hdf5-parallel-pgi/include
#xt5r4:  LIBS = -L/opt/cray/hdf5-parallel/1.8.3.0/hdf5-parallel-pgi/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm
#xt5r4:  OPT1 = -r4 -i4 \
#			  -O4 -fast -fastsse -Mipa=fast,inline -tp istanbul-32 -Minline=reshape -Mvect=nosizelimit -Munroll -Mprefetch=n:80 -Mprefetch=w -Mflushz -Minform=inform \
#			  -module $(DST) -Mpreprocess -Ktrap=divz,inv,ovf -Mdaz \
#			  -DMPI_REAL8=MPI_REAL4 -DH5T_NATIVE_DOUBLE=H5T_NATIVE_REAL \
#			  -DDGBTRF=SGBTRF -DDGBTRS=SGBTRS -DDSWAP=SSWAP -DDSCAL=SSCAL -DDGER=SGER -DDCOPY=SCOPY -DDLASWP=SLASWP -DDTRSM=STRSM -DDGEMM=SGEMM -DDCOPY=SCOPY
#xt5r4:  OPT2 = $(OPT1)
#xt5r4:  build
##-------------------------------------------------------------------------------------------------------------
#p5:     COMP = mpfort -compiler xlf90_r
#p5:     DST  = ../prog
#p5:     PRE  = -DALLOC
#p5:     CUSR = FALSE
#p5:     INCL = -I/apps/hdf5/1.6.7/parallel/lib
#p5:     LIBS = -L/apps/hdf5/1.6.7/parallel/lib -lhdf5_fortran -lhdf5 -lhdf5_hl -L/apps/zlib/1.2.3/XL/lib -lz  \
#			  /apps/szip/2.1/XL/lib/libsz.a /apps/lapack/3.1.1/XL/lib/liblapack_ppc64.a -lessl
#p5:     OPT1 = -qrealsize=8 -qintsize=4 -q64 \
#			  -O5 -qarch=pwr5 -qtune=pwr5 -qnoipa -qflttrap=zerodivide \
#			  -qmoddir=$(DST) -I$(DST) -qsuffix=cpp=f90 -WF,-Dflush=flush_
#p5:     OPT2 = $(OPT1)
#p5:     build
##-------------------------------------------------------------------------------------------------------------
##

include compile
#
