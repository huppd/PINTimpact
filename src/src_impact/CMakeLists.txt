# Include the directory itself as a path to include directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
add_definitions(-DALLOC)

# Create a variable called *_SRCS containing all .f90 files:
FILE(GLOB mod_SRCS mod_*.f90)
FILE(GLOB cmod_SRCS cmod_*.f90)
FILE(GLOB usr_SRCS usr_*.f90) 

#MESSAGE("   cmod = ${cmod}")
#MESSAGE("   cmod_SRCS = ${cmod_SRCS}")



set(limpact_srcs ${mod_SRCS} ${cmod_SRCS} ${usr_SRCS} sub_other.f90)

add_library(limpact ${limpact_srcs} )

#target_link_libraries(limpact cmod)

add_executable(impact impact.f90)

add_test( run_impact mpirun -np 4 impact)

set_target_properties(impact limpact PROPERTIES 
	COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}")

#TARGET_LINK_LIBRARIES(impact limpact hdf5 hdf5_fortran hdf5_hl z m lapack blas) 
TARGET_LINK_LIBRARIES(impact limpact ${Limpact_TPL_LIBRARIES})
