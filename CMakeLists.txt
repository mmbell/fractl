# cmake for fractl

cmake_minimum_required (VERSION 2.8)
project ( fractl )

set ( CMAKE_BUILD_TYPE Release )

# Where to find *Find*.cmake and *-config.cmake

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/" )
set ( CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH} )

find_package ( Lrose COMPONENTS Radx kd tdrp netcdf toolsa REQUIRED )	
find_package ( LibZip REQUIRED )
find_package ( GeographicLib REQUIRED )
find_package ( Eigen3 REQUIRED )

# Set RPATH so that executable finds lrose libraries without needing LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_BUILD_WITH_INSTALL_RPATH TRUE )

set ( CMAKE_INSTALL_RPATH ${LROSE_PREFIX}/lib )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# Use RUNPATH, not RPATH

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set ( CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags" )
endif()

add_definitions ( -Wall -std=c++11 )

# Specify source files

file ( GLOB sources *.cc *.hh )
set ( sources ${sources} Params.cc )	 # Add the 2 generated files

include_directories ( /usr/local/include ${LROSE_PREFIX}/include ${EIGEN3_INCLUDE_DIR} )

# Generate the argument files

link_directories( /usr/local/lib ${LROSE_PREFIX}/lib /usr/lib/x86_64-linux-gnu/hdf5/serial )
add_executable ( ${PROJECT_NAME} ${sources} )

target_link_libraries ( ${PROJECT_NAME} ${GeographicLib_LIBRARIES} ${LROSE_LIBRARIES} bz2 )

install ( TARGETS ${PROJECT_NAME} DESTINATION ${LROSE_PREFIX}/bin )
