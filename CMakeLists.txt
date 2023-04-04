cmake_minimum_required(VERSION 3.9)
project(Distance_SCAN_SIGMOD)

set(CMAKE_CXX_STANDARD 14)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(BOOST_ROOT ../boost_1_75_0)
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include)

set(Boost_NO_BOOST_CMAKE ON)

find_package(Boost REQUIRED unit_test_framework filesystem system date_time program_options serialization)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(Distance_SCAN_SIGMOD main.cpp src/config.cpp src/graph.cpp src/mylib.cpp include/query.h src/clusteralgo.cpp include/clusteralgo.h src/sketches.cpp include/sketches.h include/fhq_treap.h src/fhq_treap.cpp include/algos.h src/binary_tree.cpp include/binary_tree.h)
    target_link_libraries(Distance_SCAN_SIGMOD ${Boost_LIBRARIES})
endif ()


