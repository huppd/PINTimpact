# Include the directory itself as a path to include directories
add_subdirectory(Fields)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories("./Fields")
include_directories("./Operators")
 
# Create a variable called helloworld_SOURCES containing all .cpp files:
#set(pimpact_SOURCES pimpact.cpp)
 
# For a large number of source files you can create it in a simpler way
# using file() function:
# file(GLOB pimpact_SOURCES *.cpp)
 
# Create an executable file called helloworld from sources:
#add_executable(pimpact pimpact.cpp)

#message("cxx compiler:  ${CMAKE_CXX_COMPILER}")

#message("impactlibs:" ${ALL_f90_SRCS})
#TARGET_LINK_LIBRARIES(pimpact  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)
#TARGET_LINK_LIBRARIES(pimpact  limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)

#TARGET_LINK_LIBRARIES( pimpact  ${Pimpact_TPL_LIBRARIES} limpact ) 
#message( "tpl_libs:  ${Pimpact_TPL_LIBRARIES} ")

# Create a variable called helloworld_SOURCES containing all .cpp files:
#set(pimpact_SOURCES pimpact.cpp)
 
# For a large number of source files you can create it in a simpler way
# using file() function:
# file(GLOB pimpact_SOURCES *.cpp)
 
# Create an executable file called helloworld from sources:
add_executable(stat_stokes stat_stokes.cxx)
TARGET_LINK_LIBRARIES( stat_stokes  ${Pimpact_TPL_LIBRARIES} ) 

# Create an executable file called helloworld from sources:
add_executable( peri_stokes peri_stokes.cxx )
TARGET_LINK_LIBRARIES( peri_stokes  ${Pimpact_TPL_LIBRARIES} )

add_executable( peri_stokes2 peri_stokes2.cxx )
TARGET_LINK_LIBRARIES( peri_stokes2  ${Pimpact_TPL_LIBRARIES} ) 

#TARGET_LINK_LIBRARIES(stat_stokes  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)
#TARGET_LINK_LIBRARIES(pimpact  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)
#TARGET_LINK_LIBRARIES(pimpact  limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)


#message( "tpl_libs:  ${Pimpact_TPL_LIBRARIES} ")

