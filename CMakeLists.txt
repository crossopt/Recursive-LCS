cmake_minimum_required(VERSION 3.10) # Проверка версии CMake
project(Recursive-LCS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Werror -pedantic -pthread")

set(SOURCES
    src/monge_matrix.cpp
    src/lcs_kernel.cpp
    src/grammar_compressed.cpp
)

include_directories(inc/)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
include_directories(${Boost_INCLUDE_DIRS})

add_library(llcs STATIC ${SOURCES})
target_link_libraries(llcs ${Boost_LIBRARIES})
add_executable(lcs main.cpp)
target_link_libraries(lcs llcs)

add_executable(time_test time_test.cpp)
target_link_libraries(time_test llcs)

add_executable(dp_test dp_main.cpp)
target_link_libraries(dp_test llcs)

add_executable(lz_test lz_main.cpp)
target_link_libraries(lz_test llcs)

#GTest
enable_testing()
include(GoogleTest)
add_subdirectory(google_test)

set(TEST_SOURCES
    test/main.cpp
    test/test_grammar_compressed.cpp
    test/test_lcs_kernel.cpp
    test/test_monge_matrix.cpp
)

add_executable(lcs_test ${TEST_SOURCES})
target_link_libraries(lcs_test llcs gtest gtest_main)
add_test(NAME lcs_test COMMAND lcs_test)
