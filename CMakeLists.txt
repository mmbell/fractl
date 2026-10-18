# cmake for fractl

cmake_minimum_required(VERSION 3.20)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Default to Release build type unless specified on command line

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
endif()

# project name

project(fractl)

# Where to find *Find*.cmake and *-config.cmake

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH})

# Finding pkg-config

find_package(PkgConfig REQUIRED)

#  hdf5 and netcdf

find_package (HDF5)
find_package (NETCDF)
if (DEFINED HDF5_hdf5_LIBRARY_RELEASE)
  get_filename_component(HDF5_INSTALL_PREFIX ${HDF5_hdf5_LIBRARY_RELEASE} DIRECTORY)
endif()
message("netCDF_INSTALL_PREFIX: ${netCDF_INSTALL_PREFIX}")
message("HDF5_INSTALL_PREFIX: ${HDF5_INSTALL_PREFIX}")
message("HDF5_C_INCLUDE_DIR: ${HDF5_C_INCLUDE_DIR}")

# find LROSE
find_package(LROSE REQUIRED)

# code for eigen3 has been copied into the fractly distribution
# so we no longer need to find the package.
# find_package(Eigen3 REQUIRED)
set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/eigen3")

# enable testing for this dir and below

enable_testing()

# Set RPATH so that executable finds lrose libraries without needing LD_LIBRARY_PATH

set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Use RUNPATH, not RPATH

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
endif()

# recurse into src directory for the build

add_subdirectory(src)  

