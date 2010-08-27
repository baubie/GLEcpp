cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
project (GLETest CXX)

find_package(Boost 1.43.0 COMPONENTS iostreams filesystem)
include(${Boost_INCLUDE_DIR})

add_definitions (-Wall -W)
set (SRC src/GLE.cpp examples/basic.cpp)
add_executable (examples/basic ${SRC})
target_link_libraries(examples/basic ${Boost_LIBRARIES})
