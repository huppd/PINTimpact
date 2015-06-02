set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB belosETI_SOURCES *.cpp)

set( belosETI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Pimpact_LinearProblem.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Pimpact_LinSolverParameter.cpp PARENT_SCOPE) 
