cmake_minimum_required(VERSION 3.15)
project(CPAM VERSION 1
        DESCRIPTION "Compressed purely-functional sequences, sets, and maps in C++, with parallelism support from parlaylib."
        LANGUAGES CXX)

include(CheckCXXCompilerFlag)
include (FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

FetchContent_Declare(parlaylib
  GIT_REPOSITORY  https://github.com/cmuparlay/parlaylib.git
  GIT_TAG         master
)
FetchContent_GetProperties(parlaylib)
if(NOT parlaylib_POPULATED)
  FetchContent_Populate(parlaylib)  
  add_subdirectory(${parlaylib_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CPAM_SOURCE_DIR}/cmake")

add_library(CPAM INTERFACE)
set(CPAM_INCLUDE_DIR "${CPAM_SOURCE_DIR}/include")
target_include_directories(CPAM INTERFACE
  $<BUILD_INTERFACE:${CPAM_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_subdirectory(examples)
