# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)


include_directories("./")
include_directories("./Space")
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
include_directories("./MultiGrid/Transfers")
include_directories( "./EigenSolver" )

add_subdirectory("./Space")
add_subdirectory("./Fields")
add_subdirectory("./Operators")
add_subdirectory("./BelosPimpact")
add_subdirectory("./MultiGrid")
add_subdirectory("./MultiGrid/Transfers")
add_subdirectory("./NOXPimpact")
 

add_library( lpimpact ${Space_SOURCES} ${field_SOURCES} ${operator_SOURCES} ${belosETI_SOURCES} ${mg_SOURCES} ${trans_SOURCES} ${nox_SOURCES} )
set_target_properties( lpimpact PROPERTIES LINKER_LANGUAGE CXX )
target_link_libraries( lpimpact cmod )
 

#message("cxx compiler:  ${CMAKE_CXX_COMPILER}")
#message("impactlibs:" ${ALL_f90_SRCS})

#add_executable(stat_stokes stat_stokes.cxx)
#TARGET_LINK_LIBRARIES( stat_stokes  ${Pimpact_TPL_LIBRARIES} ) 

#add_executable( peri_stokes peri_stokes.cxx )
#TARGET_LINK_LIBRARIES( peri_stokes  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff Pimpact_ConvDiff.cxx )
TARGET_LINK_LIBRARIES( conv_diff ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff_mg Pimpact_ConvDiffMG.cxx )
TARGET_LINK_LIBRARIES( conv_diff_mg  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff_mg2 Pimpact_ConvDiffMG2.cxx )
TARGET_LINK_LIBRARIES( conv_diff_mg2  ${Pimpact_TPL_LIBRARIES} )

add_executable( conv_diff2 Pimpact_ConvDiff2.cxx )
TARGET_LINK_LIBRARIES( conv_diff2  ${Pimpact_TPL_LIBRARIES} )

add_executable( peri_navier3D peri_navier3D.cxx )
TARGET_LINK_LIBRARIES( peri_navier3D ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navierXML peri_navierSMGXML.cxx )
TARGET_LINK_LIBRARIES( peri_navierXML ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier4 peri_navier4.cxx )
TARGET_LINK_LIBRARIES( peri_navier4 ${Pimpact_TPL_LIBRARIES} ) 
