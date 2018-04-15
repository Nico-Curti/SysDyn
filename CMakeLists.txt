cmake_minimum_required (VERSION 3.10)
project (SystemDynamics)
set (CMAKE_CXX_STANDARD 14)

# make sure that the default is a RELEASE
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

### Set output directories on a per-configuration base
# Single configuration
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )
# Multi configuration
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_SOURCE_DIR}/bin )
endforeach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )

###########################################################################################
add_executable(ChemicalMasterEquation   ${CMAKE_SOURCE_DIR}/cpp/ChemicalMasterEquation.cpp)
add_executable(diffusion2D              ${CMAKE_SOURCE_DIR}/cpp/diffusion2D.cpp)
#add_executable(FalseNearesetNeighbors   ${CMAKE_SOURCE_DIR}/cpp/FalseNearesetNeighbors.cpp)
add_executable(GaussSeidelSolve         ${CMAKE_SOURCE_DIR}/cpp/GaussSeidelSolve.cpp)
add_executable(GrassbergerProcaccia     ${CMAKE_SOURCE_DIR}/cpp/GrassbergProcaccia.cpp)
add_executable(SpringLayout             ${CMAKE_SOURCE_DIR}/cpp/SpringLayout.cpp)
add_executable(ThomasSolve              ${CMAKE_SOURCE_DIR}/cpp/ThomasSolve.cpp)

