# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
#add_definitions(-DALLOC -xf95-cpp-input)
add_definitions(-DALLOC)

# Create a variable called helloworld_SOURCES containing all .cpp files:
#set(pimpact_SOURCES *.f90)
FILE(GLOB mod_SRCS mod_*.f90)
FILE(GLOB cmod_SRCS cmod_*.f90)
FILE(GLOB usr_SRCS usr_*.f90) 
set(limpact_srcs ${mod_SRCS} ${cmod_SRCS} ${usr_SRCS} sub_other.f90)

# For a large number of source files you can create it in a simpler way
# using file() function:
# file(GLOB pimpact_SOURCES *.cpp)

#MESSAGE("Fortran compiler: " ${CMAKE_Fortran_COMPILER})
#Message("Fortran_FORMAT: " ${Fortran_FORMAT})
#MESSAGE("CMAKE_Fortran_COMPILE_FLAGS: " ${CMAKE_Fortran_COMPILE_FLAGS})
#Message("Impact libs: " ${ALL_f90_SRCS})

add_library(limpact ${limpact_srcs} )
#add_library(mod ${mod_SRCS})
#add_library(usr ${usr_SRCS})

add_executable(impact impact.f90)

add_test( run_impact mpirun -np 4 impact)

set_target_properties(impact limpact PROPERTIES 
	COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}")

#TARGET_LINK_LIBRARIES(impact limpact hdf5 hdf5_fortran hdf5_hl z m lapack blas) 
TARGET_LINK_LIBRARIES(impact limpact ${Limpact_TPL_LIBRARIES}) 
