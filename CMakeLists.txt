cmake_minimum_required (VERSION 3.8.2)
project (SysDyn LANGUAGES CXX VERSION 1.0.0 DESCRIPTION "System Dynamics Functions and Examples")
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# SysDyn Version
set (MAJOR    1)
set (MINOR    0)
set (REVISION 0)
add_definitions (-DMAJOR=${MAJOR} -DMINOR=${MINOR} -DREVISION=${REVISION})

#################################################################
#                         COMPILE OPTIONS                       #
#################################################################


#################################################################
#                         SETTING VARIABLES                     #
#################################################################

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})

if ( NOT APPLE )
  set (CMAKE_SKIP_BUILD_RPATH             FALSE )
  set (CMAKE_BUILD_WITH_INSTALL_RPATH     FALSE )
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH  TRUE  )
endif()

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
#                         PARSE OPTIONS                         #
#################################################################

if (FORCE_LOCAL_OPENCV)
  #find_package (OpenCV COMPONENTS video text core imgproc imgcodecs highgui plot REQUIRED PATHS ${OpenCV_PATH} NO_DEFAULT_PATH)
  find_package (OpenCV REQUIRED PATHS ${OpenCV_PATH} NO_DEFAULT_PATH)
else ()
  if (APPLE)
    #find_package (OpenCV COMPONENTS video text core imgproc imgcodecs highgui plot REQUIRED PATHS /usr/local/opt NO_DEFAULT_PATH)
    find_package (OpenCV REQUIRED PATHS /usr/local/opt NO_DEFAULT_PATH)
  else ()
    #find_package (OpenCV COMPONENTS video text core imgproc imgcodecs highgui plot REQUIRED)
    find_package (OpenCV REQUIRED)
  endif ()
endif ()

if (OpenCV_FOUND)
  message (STATUS "OpenCV found: ${OpenCV_INCLUDE_DIRS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENCV")
  include_directories (${OpenCV_INCLUDE_DIRS})

  list(APPEND linked_libs ${OpenCV_LIBS})
endif ()

message (STATUS "Looking for pthread library")
set (CMAKE_THREAD_PREFER_PTHREAD ON)
find_package (Threads REQUIRED)
list (APPEND linked_libs Threads::Threads )


if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
endif()


set (OUT_DIR      ${CMAKE_SOURCE_DIR}/bin             CACHE PATH "Path where outputs will be installed"        FORCE)

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

add_executable(ChemicalMasterEquation   ${CMAKE_SOURCE_DIR}/cpp/ChemicalMasterEquation.cpp)
target_link_libraries(ChemicalMasterEquation ${linked_libs})

add_executable(brusselator              ${CMAKE_SOURCE_DIR}/cpp/brusselator.cpp)
target_link_libraries(brusselator ${linked_libs})

add_executable(bernoulli2D              ${CMAKE_SOURCE_DIR}/cpp/bernoulli2D.cpp)
target_link_libraries(bernoulli2D ${linked_libs})

add_executable(GaussSeidelSolve         ${CMAKE_SOURCE_DIR}/cpp/GaussSeidelSolve.cpp)
target_link_libraries(GaussSeidelSolve ${linked_libs})

add_executable(GrassbergerProcaccia     ${CMAKE_SOURCE_DIR}/cpp/GrassbergProcaccia.cpp)

add_executable(SpringLayout             ${CMAKE_SOURCE_DIR}/cpp/SpringLayout.cpp)
target_link_libraries(SpringLayout ${linked_libs})

add_executable(ThomasSolve              ${CMAKE_SOURCE_DIR}/cpp/ThomasSolve.cpp)
target_link_libraries(ThomasSolve ${linked_libs})

#################################################################
#                          INSTALLERS                           #
#################################################################

install(TARGETS ChemicalMasterEquation      DESTINATION ${OUT_DIR})
install(TARGETS brusselator                 DESTINATION ${OUT_DIR})
install(TARGETS bernoulli2D                 DESTINATION ${OUT_DIR})
install(TARGETS GaussSeidelSolve            DESTINATION ${OUT_DIR})
install(TARGETS GrassbergerProcaccia        DESTINATION ${OUT_DIR})
install(TARGETS SpringLayout                DESTINATION ${OUT_DIR})
install(TARGETS ThomasSolve                 DESTINATION ${OUT_DIR})