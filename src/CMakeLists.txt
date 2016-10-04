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

find_package( Trilinos REQUIRED )


# Echo trilinos build info just for fun
message( "\nFound Trilinos!  Here are the details: " )
message( "Trilinos_DIR = ${Trilinos_DIR}" )
message( "Trilinos_VERSION = ${Trilinos_VERSION}" )
message( "Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}" )
message( "Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}" )
message( "Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}" )
message( "Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}" )
message( "Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}" )
message( "Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}" )
message( "Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}" )
message( "Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}" )
message( "Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}" )
message( "End of Trilinos details\n" )



find_package( MPI REQUIRED )
message( "\nMPI CXX? ${MPI_CXX_FOUND}" )
message( "MPI CXX INCLUDE_DIRS = ${MPI_CXX_INCLUDE_PATH}" )
message( "MPI CXX LIBRARIES =    ${MPI_CXX_LIBRARIES}" )
message( "\nMPI Fortran? ${MPI_Fortran_FOUND}" )
message( "MPI Fortran INCLUDE_DIRS = ${MPI_Fortran_INCLUDE_PATH}" )
message( "MPI Fortran LIBRARIES =    ${MPI_Fortran_LIBRARIES}\n" )

option( findHDF5 "findHDF5" OFF)

if( findHDF5 )
	FIND_PACKAGE( HDF5 COMPONENTS C Fortran HL Fortran_HL )
	message( "\nHDF5? ${HDF5_FOUND}" )
	message( "HDF5_IS_PARALLEL? ${HDF5_IS_PARALLEL}" )
	message( "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}" )
	message( "HDF5_DEFINITIONS: ${HDF5_DEFINITIONS}" )
	message( "HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}" )
	message( "HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}" )
	message( "HDF5_Fortran_LIBRARIES: ${HDF5_Fortran_LIBRARIES}" )
	message( "HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}" )
	message( "HDF5_Fortran_HL_LIBRARIES: ${HDF5_Fortran_HL_LIBRARIES}" )
	message( "HDF5_LIBRARIES: ${HDF5_LIBRARIES}" )
	message( "HDF5_LIBRARY_DIRS: ${HDF5_LIBRARY_DIRS}" )

	add_definitions( ${HDF5_DEFINITIONS} )
else()
	message( "findHDF5? off ")
endif()


option( enable_ETI "enable_ETI" OFF )

if( enable_ETI )
	add_definitions( -DCOMPILE_ETI )
	message( "Pimpact:	enabled ETI\n" )
else()
	message( "Pimpact:	disabled ETI\n" )
endif()

# Make sure to use same compilers and flags as Trilinos
set( CMAKE_CXX_COMPILER     ${Trilinos_CXX_COMPILER}     )
set( CMAKE_C_COMPILER       ${Trilinos_C_COMPILER}       )
set( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )


set( CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS}     ${Trilinos_CXX_COMPILER_FLAGS}     -Wall -pedantic -std=c++11" )
set( CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${Trilinos_C_COMPILER_FLAGS}       -Wall -pedantic" )
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -Wall" )


if( CMAKE_COMPILER_IS_GNUCXX )
	set( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -ftree-vectorizer-verbose=0 -funroll-all-loops -ftree-vectorize  -ffast-math -funsafe-math-optimizations -fargument-noalias-global -march=native -malign-double ")
	set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} ${CMAKE_CXX_FLAGS_RELEASE}  -g ")# -pg
	set( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fdollar-ok -Wall -cpp " )
	set( CMAKE_Fortran_FLAGS_RELEASE  " ${CMAKE_Fortran_FLAGS_RELEASE} -Ofast -ftree-vectorizer-verbose=0  -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno  -Ofast -fargument-noalias-global -march=native -malign-double" )
	set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} ${CMAKE_Fortran_FLAGS_RELEASE}  -g  ")  # -pg
	set( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -xf95-cpp-input -fdefault-real-8 -fdefault-double-8  "  )
else()
	set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -vec-report2 -ip -fno-alias -march=native " )
	set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} ${CMAKE_CXX_FLAGS_RELEASE}  -g ")# -pg 
	set( CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_Fortran_FLAGS_RELEASE} -funroll-all-loops -march=native -fomit-frame-pointer -fno-math-errno -vec-report2 -ip -fno-alias  " )
	set( CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} -free -fpp -r8   "  )
	set( CMAKE_Fortran_COMPILE_FLAGS_RELEASE "${CMAKE_Fortran_COMPILE_FLAGS_RELEASE} ${Trilinos_Fortran_COMPILER_FLAGS} -ip -fno-alias -march=native "  )
enDIF()


message( "Fortran compiler: ${CMAKE_Fortran_COMPILER}" )
message( "Fortran flags: " ${CMAKE_Fortran_FLAGS})
message( "\nCXX compiler: ${CMAKE_CXX_COMPILER}" )
message( "CXX flags: " ${CMAKE_CXX_FLAGS})
message( "CMAKE_Fortran_COMPILE_FLAGS: ${CMAKE_Fortran_COMPILE_FLAGS}\n" )


#if( findHDF5 )
	#set( Limpact_TPL_LIBRARIES ${HDF5_LIBRARIES} ${Trilinos_TPL_LIBRARIES} )
#else()
	set( Limpact_TPL_LIBRARIES hdf5hl_fortran hdf5_hl hdf5_fortran hdf5 z m
		${Trilinos_TPL_LIBRARIES} )
#endif()

set( Pimpact_TPL_LIBRARIES lpimpact ${Limpact_TPL_LIBRARIES} gfortran
	${MPI_Fortran_LIBRARIES}  ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} )

message( "\nLimpact_TPL_LIBRARIES: ${Limpact_TPL_LIBRARIES}\n" )
message( "\nPimpact_TPL_LIBRARIES: ${Pimpact_TPL_LIBRARIES}\n" )

include_directories( ${inc_dir} ${Trilinos_INCLUDE_DIRS}
	${Trilinos_TPL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
link_directories( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
	${HDF5_LIBRARY_DIRS} )


# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (ie. MacroOutOfSourceBuild.cmake):
set( CMAKE_MODULE_PATH ${pimpact_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH} )

get_property( inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES )
message( "Include dirs = ${inc_dirs}" )
get_property( link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES )
message( "Link dirs = ${link_dirs}\n" )

enable_testing()



add_subdirectory( src_impact )
add_subdirectory( src_f )
add_subdirectory( src_c )
add_subdirectory( test )
