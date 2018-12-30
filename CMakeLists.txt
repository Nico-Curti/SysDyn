cmake_minimum_required (VERSION 3.8.2)
project (SysDyn LANGUAGES CXX VERSION 1.0.0 DESCRIPTION "System Dynamics Functions and Examples")
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# SysDyn Version
set (MAJOR    1)
set (MINOR    0)
set (REVISION 0)

if ( CMAKE_COMPILER_IS_GNUCC )
  add_compile_options (-Wall -Wextra -Wno-unused-result)
  string (REGEX REPLACE "-O3" "-Ofast" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
endif()
if ( MSVC )
  add_compile_options (/wd4028)
  add_compile_options (/wd4244)
  add_compile_options (/wd4267)
  add_compile_options (/wd4305)
  add_compile_options (/wd4477)
  add_compile_options (/wd4996)
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRTD")
  #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRT")
  set (CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

#################################################################
#                          SUMMARY                              #
#################################################################

message(STATUS ""                                                                     )
message(STATUS "=================== SysDyn configuration Summary =================="  )
message(STATUS "   SysDyn version: ${MAJOR}.${MINOR}.${REVISION}"                     )
message(STATUS ""                                                                     )
message(STATUS "   C++ :"                                                             )
message(STATUS "      C++ Compiler : ${CMAKE_CXX_COMPILER}"                           )
message(STATUS "      C++ flags    :"                                                 )
foreach(FLAG ${CMAKE_CXX_FLAGS})
  message(STATUS "                    * ${FLAG}"                                      )
endforeach(FLAG)
foreach(FLAG ${CMAKE_CXX_FLAGS_RELEASE})
  message(STATUS "                    * ${FLAG}"                                      )
endforeach(FLAG)
message(STATUS "      Linker flags : "                                                )
foreach(FLAG ${linked_libs})
  message(STATUS "                    * ${FLAG}"                                      )
endforeach(FLAG)
message(STATUS ""                                                                     )

#################################################################
#                         MAIN RULES                            #
#################################################################

add_custom_target(ChemicalMasterEquation_)
add_executable(ChemicalMasterEquation   ${CMAKE_SOURCE_DIR}/cpp/ChemicalMasterEquation.cpp)

add_custom_target(diffusion2D_)
add_executable(diffusion2D              ${CMAKE_SOURCE_DIR}/cpp/diffusion2D.cpp)

add_custom_target(FalseNearestNeighbors_)
add_executable(FalseNearestNeighbors   ${CMAKE_SOURCE_DIR}/cpp/FalseNearestNeighbors.cpp)

add_custom_target(GaussSeidelSolve_)
add_executable(GaussSeidelSolve         ${CMAKE_SOURCE_DIR}/cpp/GaussSeidelSolve.cpp)

add_custom_target(GrassbergerProcaccia_)
add_executable(GrassbergerProcaccia     ${CMAKE_SOURCE_DIR}/cpp/GrassbergProcaccia.cpp)

add_custom_target(SpringLayout_)
add_executable(SpringLayout             ${CMAKE_SOURCE_DIR}/cpp/SpringLayout.cpp)

add_custom_target(ThomasSolve_)
add_executable(ThomasSolve              ${CMAKE_SOURCE_DIR}/cpp/ThomasSolve.cpp)

#################################################################
#                          INSTALLERS                           #
#################################################################

install(TARGETS ChemicalMasterEquation      DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS diffusion2D                 DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS FalseNearestNeighbors       DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS GaussSeidelSolve            DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS GrassbergerProcaccia        DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS SpringLayout                DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})
install(TARGETS ThomasSolve                 DESTINATION ${CMAKE_SOURCE_DIR}/${OUT_DIR})