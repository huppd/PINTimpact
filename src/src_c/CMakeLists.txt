# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)


include_directories("./")
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
include_directories("./MultiGrid/Transfers")
include_directories( "./EigenSolver" )

add_subdirectory("./Grid")
add_subdirectory("./Fields")
add_subdirectory("./Operators")
add_subdirectory("./BelosPimpact")
add_subdirectory("./MultiGrid")
add_subdirectory("./MultiGrid/Transfers")
add_subdirectory("./NOXPimpact")
 

add_library( lpimpact ${Grid_SOURCES} ${field_SOURCES} ${operator_SOURCES} ${belosETI_SOURCES} ${mg_SOURCES} ${trans_SOURCES} ${nox_SOURCES} )
set_target_properties( lpimpact PROPERTIES LINKER_LANGUAGE CXX )
target_link_libraries( lpimpact cmod )
 

#message("cxx compiler:  ${CMAKE_CXX_COMPILER}")
#message("impactlibs:" ${ALL_f90_SRCS})

add_executable( convDiff convDiff.cxx )
TARGET_LINK_LIBRARIES( convDiff ${Pimpact_TPL_LIBRARIES} )

add_executable( convDiff_mg convDiffMG.cxx )
TARGET_LINK_LIBRARIES( convDiff_mg  ${Pimpact_TPL_LIBRARIES} )

add_executable( convDiff_mg2 convDiffMG2.cxx )
TARGET_LINK_LIBRARIES( convDiff_mg2  ${Pimpact_TPL_LIBRARIES} )

add_executable( convDiff2 convDiff2.cxx )
TARGET_LINK_LIBRARIES( convDiff2  ${Pimpact_TPL_LIBRARIES} )

add_executable( modeConvDiff modeConvDiff.cxx )
TARGET_LINK_LIBRARIES( modeConvDiff ${Pimpact_TPL_LIBRARIES} )

add_executable( peri_navier2D peri_navier2D.cxx )
TARGET_LINK_LIBRARIES( peri_navier2D ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier3D peri_navier3D.cxx )
TARGET_LINK_LIBRARIES( peri_navier3D ${Pimpact_TPL_LIBRARIES} ) 

add_executable( analyzer3D analyzer3D.cxx )
TARGET_LINK_LIBRARIES( analyzer3D ${Pimpact_TPL_LIBRARIES} ) 

add_executable( peri_navier3DTime peri_navier3DTime.cxx )
TARGET_LINK_LIBRARIES( peri_navier3DTime ${Pimpact_TPL_LIBRARIES} ) 
