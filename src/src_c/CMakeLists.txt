# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories("./Space")
include_directories("./Grid")
include_directories("./Fields")
include_directories("./Operators")
include_directories("./Operators/BasicOperators")
include_directories("./Operators/CompoundOperators")
include_directories("./Operators/ModeOperators")
include_directories("./Operators/MultiHarmonicOperators")
include_directories("./Operators/TimeOperators")
include_directories("./BelosPimpact")
include_directories("./NOXPimpact")
include_directories("./MultiGrid")

 
 
#message("cxx compiler:  ${CMAKE_CXX_COMPILER}")
#message("impactlibs:" ${ALL_f90_SRCS})

# Create an executable file called helloworld from sources:
add_executable(stat_stokes stat_stokes.cxx)
TARGET_LINK_LIBRARIES( stat_stokes  ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_stokes peri_stokes.cxx )
TARGET_LINK_LIBRARIES( peri_stokes  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff Pimpact_ConvDiff.cxx )
TARGET_LINK_LIBRARIES( conv_diff  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff_mg Pimpact_ConvDiffMG.cxx )
TARGET_LINK_LIBRARIES( conv_diff_mg  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff2 Pimpact_ConvDiff2.cxx )
TARGET_LINK_LIBRARIES( conv_diff2  ${Pimpact_TPL_LIBRARIES} )

#add_executable( peri_navier_fix peri_navier_fix.cxx )
#TARGET_LINK_LIBRARIES( peri_navier_fix ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier peri_navier.cxx )
TARGET_LINK_LIBRARIES( peri_navier ${Pimpact_TPL_LIBRARIES} ) 

#add_executable( peri_navier4 peri_navier4.cxx )
#TARGET_LINK_LIBRARIES( peri_navier4 ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier_timing peri_navier_timing.cxx )
TARGET_LINK_LIBRARIES( peri_navier_timing ${Pimpact_TPL_LIBRARIES} ) 

