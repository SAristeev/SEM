cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)
project(SEM VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include_directories("external")

find_package(VTK REQUIRED)

if(VTK_FOUND)
  message("-- Found VTK. Version:" ${VTK_VERSION})
  message("-- VTK_DIR: ${VTK_DIR}")
else(VTK_FOUND)
  message(FATAL_ERROR "Cannot build the executable without VTK. Please set the VTK_DIR")
endif(VTK_FOUND)

find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
find_package(TBB CONFIG REQUIRED PATHS $ENV{MKLROOT})

add_subdirectory(src)