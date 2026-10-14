cmake_minimum_required(VERSION 3.15)
project(svaha VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math")
    # OpenMP on macOS with Clang often requires extra work, 
    # but we are using C++17 threads/async which don't need -fopenmp.
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -fopenmp")
endif()

# Include directories
include_directories(include)

add_executable(svaha src/main.cpp)

# Target properties
target_include_directories(svaha PUBLIC include)

# Tests
add_executable(test_svaha tests/test_svaha.cpp)
target_include_directories(test_svaha PUBLIC include)
enable_testing()
add_test(NAME svaha_unit_tests COMMAND test_svaha)
