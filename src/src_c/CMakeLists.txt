# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories("./Space")
include_directories("./Grid")
include_directories("./Fields")
include_directories("./Operators")
include_directories("./NOX_Pimpact")
include_directories("./MultiGrid")

 
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

add_executable( peri_stokes peri_stokes.cxx )
TARGET_LINK_LIBRARIES( peri_stokes  ${Pimpact_TPL_LIBRARIES} )

#add_executable( peri_stokes2 peri_stokes2.cxx )
#TARGET_LINK_LIBRARIES( peri_stokes2  ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_burgers peri_burgers.cxx )
TARGET_LINK_LIBRARIES( peri_burgers  ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_burgers_timing peri_burgers_timing.cxx )
TARGET_LINK_LIBRARIES( peri_burgers_timing  ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier peri_navier.cxx )
TARGET_LINK_LIBRARIES( peri_navier ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier4 peri_navier4.cxx )
TARGET_LINK_LIBRARIES( peri_navier4 ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier_timing peri_navier_timing.cxx )
TARGET_LINK_LIBRARIES( peri_navier_timing ${Pimpact_TPL_LIBRARIES} ) 

#add_executable( mlguide mlguide.cxx )
#TARGET_LINK_LIBRARIES( mlguide ${Pimpact_TPL_LIBRARIES} ) 

#add_executable( mlguide_par mlguide_par.cxx )
#TARGET_LINK_LIBRARIES( mlguide_par ${Pimpact_TPL_LIBRARIES} ) 

#TARGET_LINK_LIBRARIES(stat_stokes  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)
#TARGET_LINK_LIBRARIES(pimpact  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)
#TARGET_LINK_LIBRARIES(pimpact  limpact hdf5 hdf5_fortran hdf5_hl z mpi_f90 mpi_f77 gfortran lapack blas m)


#message( "tpl_libs:  ${Pimpact_TPL_LIBRARIES} ")
