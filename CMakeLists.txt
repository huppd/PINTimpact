# CMAKE File for "MyApp" application building against an installed Trilinos

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

cmake_minimum_required(VERSION 2.8)

# Use Trilinos_PREFIX, if the user set it, to help find Trilinos.
# The final location will actually be held in Trilinos_DIR which must
# point at "<prefix>/lib/cmake/Trilinos", but this helps the search.
#SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
#SET(CMAKE_PREFIX_PATH "/usr/lib/")
#message("CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH})

# Get Trilinos as one entity
FIND_PACKAGE(Trilinos REQUIRED)
#FIND_PACKAGE(HDF5 NAMES hdf5 REQUIRED)
#FIND_PACKAGE(HDF5 REQUIRED)


# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

PROJECT(pimpact CXX Fortran)

# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )


SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS} -std=c++0x -Wall -pedantic")
SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}" )

#if(  "${CMAKE_BUILD_TYPE}" equal "RELEASE" )
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -funroll-all-loops -ftree-vectorize -march=native ")
	SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS}  -ffree-line-length-none  -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno -fdollar-ok -Wall" )
	SET( CMAKE_Fortran_COMPILE_FLAGS "${Trilinos_Fortran_COMPILER_FLAGS} -xf95-cpp-input -fdefault-real-8 -fdefault-double-8 "  )
ELSE()
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST " )
	SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -funroll-all-loops -march=native -fomit-frame-pointer -fno-math-errno -xHOST " )
	SET( CMAKE_Fortran_COMPILE_FLAGS "${Trilinos_Fortran_COMPILER_FLAGS} -free -fpp -r8 "  )
	#SET( CMAKE_Fortran_LINK_FLAGS "${Trilinos_Fortran_LINK_FLAGS} -static "  )
ENDIF()
#ENDIF()

MESSAGE("Fortran flags: " ${CMAKE_Fortran_FLAGS})
MESSAGE("CXX flags: " ${CMAKE_CXX_FLAGS})


set( Limpact_TPL_LIBRARIES hdf5hl_fortran hdf5_hl hdf5_fortran hdf5  z m lapack blas )
set( Pimpact_TPL_LIBRARIES limpact ${Limpact_TPL_LIBRARIES} gfortran mpi_f77 mpi_f90   ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}   )

MESSAGE("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
MESSAGE("CMAKE_Fortran_COMPILE_FLAGS:  ${CMAKE_Fortran_COMPILE_FLAGS}")

#ADD_DEFINITIONS(-DMYAPP_EPETRA)

INCLUDE_DIRECTORIES(${inc_dir} ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

#find_path(STDIO_INCLUDE_PATH stdio.h)
#include_directories("${STDIO_INCLUDE_PATH}/dummy/../")

#ADD_LIBRARY(myappLib src_file.cpp src_file.hpp)

#ADD_EXECUTABLE(MyApp.exe main_file.cpp)

#TARGET_LINK_LIBRARIES(MyApp.exe  myappLib  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 

#enable_testing()
#add_test(NAME MyTest COMMAND MyApp.exe)

# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable which stores the
# directories of additional CMake modules (ie. MacroOutOfSourceBuild.cmake):
set(CMAKE_MODULE_PATH ${pimpact_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})


get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("Include dirs = ${inc_dirs}")
get_property(link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES)
message("Link dirs = ${link_dirs}")

enable_testing()


add_subdirectory(src)
add_subdirectory(test)
