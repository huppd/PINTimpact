set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB Space_SOURCES *.cpp)

add_library( Space ${Space_SOURCES} )
#add_dependencies( Space cmod )
target_link_libraries(Space cmod)
