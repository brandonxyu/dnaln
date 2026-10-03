# Optional CMake build (the Makefile is the primary, tested build).
#   cmake -B build-cmake -DCMAKE_BUILD_TYPE=Release && cmake --build build-cmake -j
#   ctest --test-dir build-cmake
cmake_minimum_required(VERSION 3.16)
project(dnaln CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(ZLIB REQUIRED)
include(CheckCXXCompilerFlag)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
  check_cxx_compiler_flag("-march=native" HAVE_MARCH_NATIVE)
  if(HAVE_MARCH_NATIVE)
    add_compile_options(-march=native)
  endif()
else()
  check_cxx_compiler_flag("-mcpu=native" HAVE_MCPU_NATIVE)
  if(HAVE_MCPU_NATIVE)
    add_compile_options(-mcpu=native)
  endif()
endif()
add_compile_options(-Wall -Wextra -Wno-unused-parameter)

add_library(dnaln_core STATIC
  src/common.cpp src/seqio.cpp src/reference.cpp src/minimizer.cpp src/index.cpp
  src/chain.cpp src/align.cpp src/align_simd.cpp src/mapper.cpp src/simulate.cpp
  src/evaluate.cpp src/groundtruth.cpp src/bench.cpp src/demo.cpp)
target_include_directories(dnaln_core PUBLIC src)
target_link_libraries(dnaln_core PUBLIC ZLIB::ZLIB)

add_executable(dnaln src/main.cpp)
target_link_libraries(dnaln PRIVATE dnaln_core)

add_executable(dnaln_tests tests/test_main.cpp)
target_link_libraries(dnaln_tests PRIVATE dnaln_core)

enable_testing()
add_test(NAME unit_tests COMMAND dnaln_tests)
