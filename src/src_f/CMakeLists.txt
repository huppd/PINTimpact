#add_definitions(-DALLOC)

# Create a variable called *_SRCS containing all .f90 files:
FILE(GLOB cmod_SRCS cmod_*.f90)

add_library( cmod ${cmod_SRCS} )

set_target_properties( cmod PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS}")
