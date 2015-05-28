set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB field_SOURCES *.cpp)

add_library( field ${field_SOURCES} )
#add_dependencies( field cmod )
target_link_libraries( field Space )
