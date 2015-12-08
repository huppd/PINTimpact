# CMAKE File for "PIMPACT" application building against an installed Trilinos

#This file was created by modifiying the files in
#Trilinos/demos/buildAgaintsTrilinos. The primary change was to make it a single
#file cmake system in a flat directory. If you would like to run a cmake
#configure using this file you should grab this file and src_file.cpp,
#src_file.hpp, main_file.cpp from buildAgainstTrilinos and place them in a new
#directory. From there you can run:
#"cmake -DTrilinos_PREFIX=<path to trilinos>." to configure. Another
#important change is the buildAgainstTrilinos does some checking to see which
#packages and tpls are enabled and behaves accordingly. However, this file does
#only a serial configure(no mpi) and assumes that the install of Trilinos it is
#pointed to has Epetra enabled. 

cmake_minimum_required(VERSION 2.8.10)

project( pimpact CXX Fortran )

# Get Trilinos as one entity
# Use Trilinos_PREFIX, if the user set it, to help find Trilinos.
# The final location will actually be held in Trilinos_DIR which must
# point at "<prefix>/lib/cmake/Trilinos", but this helps the search.

find_package(Trilinos REQUIRED)


# Echo trilinos build info just for fun
message( "\nFound Trilinos!  Here are the details: " )
message( "   Trilinos_DIR = ${Trilinos_DIR}" )
message( "   Trilinos_VERSION = ${Trilinos_VERSION}" )
message( "   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}" )
message( "   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}" )
message( "   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}" )
message( "   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}" )
message( "   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}" )
message( "   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}" )
message( "   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}" )
message( "   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}" )
message( "   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}" )
message( "End of Trilinos details\n" )




find_package( MPI REQUIRED )
message( "/nMPI CXX? ${MPI_CXX_FOUND}" )
message( "MPI CXX INCLUDE_DIRS = ${MPI_CXX_INCLUDE_PATH}" )
message( "MPI CXX LIBRARIES =    ${MPI_CXX_LIBRARIES}" )
message( "\nMPI Fortran? ${MPI_Fortran_FOUND}" )
message( "MPI Fortran INCLUDE_DIRS = ${MPI_Fortran_INCLUDE_PATH}" )
message( "MPI Fortran LIBRARIES =    ${MPI_Fortran_LIBRARIES}\n" )

#FIND_PACKAGE(HDF5 NAMES hdf5 REQUIRED)
#FIND_PACKAGE(HDF5 REQUIRED)

option( enable_ETI "enable_ETI" OFF)

if( enable_ETI )
	add_definitions( -DCOMPILE_ETI )
	message( "Pimpact:	enabled ETI" )
else()
	message( "Pimpact:	disabled ETI" )
endif()

# Make sure to use same compilers and flags as Trilinos
set( CMAKE_CXX_COMPILER     ${Trilinos_CXX_COMPILER}     )
set( CMAKE_C_COMPILER       ${Trilinos_C_COMPILER}       )
set( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )


set( CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS}     ${Trilinos_CXX_COMPILER_FLAGS}     -Wall -pedantic -std=c++11" )
set( CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${Trilinos_C_COMPILER_FLAGS}       -Wall -pedantic" )
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -Wall" )


IF( CMAKE_COMPILER_IS_GNUCXX )
	SET( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -ftree-vectorizer-verbose=1 -funroll-all-loops -ftree-vectorize  -ffast-math -funsafe-math-optimizations -fargument-noalias-global -march=native -malign-double ")
	set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} ${CMAKE_CXX_FLAGS_RELEASE}  -g ")# -pg
	SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fdollar-ok -Wall -cpp " )
	SET( CMAKE_Fortran_FLAGS_RELEASE  " ${CMAKE_Fortran_FLAGS_RELEASE} -ftree-vectorizer-verbose=1  -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno  -Ofast -fargument-noalias-global -march=native -malign-double" )
	set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} ${CMAKE_Fortran_FLAGS_RELEASE}  -g  ")  # -pg
	SET( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -xf95-cpp-input -fdefault-real-8 -fdefault-double-8  "  )
ELSE()
	SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -vec-report2 -ip -fno-alias -march=native " )
	SET( CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_Fortran_FLAGS_RELEASE} -funroll-all-loops -march=native -fomit-frame-pointer -fno-math-errno -vec-report2 -ip -fno-alias  " )
	SET( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -free -fpp -r8   "  )
	SET( CMAKE_Fortran_COMPILE_FLAGS_RELEASE "${CMAKE_Fortran_COMPILE_FLAGS_RELEASE} ${Trilinos_Fortran_COMPILER_FLAGS} -ip -fno-alias -march=native "  )
	set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} ${CMAKE_CXX_FLAGS_RELEASE}  -g ")# -pg 
ENDIF()


message("Fortran flags: " ${CMAKE_Fortran_FLAGS})
message("CXX flags: " ${CMAKE_CXX_FLAGS})


set( Limpact_TPL_LIBRARIES hdf5hl_fortran hdf5_hl hdf5_fortran hdf5 z m ${Trilinos_TPL_LIBRARIES} )
#set( Pimpact_TPL_LIBRARIES lpimpact ${Limpact_TPL_LIBRARIES} gfortran mpi_f77 mpi_f90  ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} )
set( Pimpact_TPL_LIBRARIES lpimpact ${Limpact_TPL_LIBRARIES} gfortran ${MPI_Fortran_LIBRARIES}  ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} )


message( "Fortran compiler: ${CMAKE_Fortran_COMPILER}" )
message( "CMAKE_Fortran_COMPILE_FLAGS: ${CMAKE_Fortran_COMPILE_FLAGS}" )


INCLUDE_DIRECTORIES( ${inc_dir} ${Trilinos_INCLUDE_DIRS}
	${Trilinos_TPL_INCLUDE_DIRS}
	"/cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/include" )
LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )


message( "Include: ${inc_dir} ${Trilinos_INCLUDE_DIRS}
${Trilinos_TPL_INCLUDE_DIRS} ${HDF5_LIBRARY_DIRS} " )
#message( "Link: ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}" )

# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (ie. MacroOutOfSourceBuild.cmake):
set( CMAKE_MODULE_PATH ${pimpact_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH} )

get_property( inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES )
message( "Include dirs = ${inc_dirs}" )
get_property( link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES )
message( "Link dirs = ${link_dirs}" )

enable_testing()



add_subdirectory( src_impact )
add_subdirectory( src_f )
add_subdirectory( src_c )
add_subdirectory( test )
