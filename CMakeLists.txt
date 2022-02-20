# This is the CMake script for compiling a GTVimageVect project.
project( GTVimageVect ) 


cmake_minimum_required (VERSION 3.22) 
cmake_policy(SET CMP0057 NEW) # the new interpretation of IN_LIST is by default not set (at least until  3.21.1) else findiing boost dependancie fails.

include_directories ("${PROJECT_SOURCE_DIR}/src")

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fext-numeric-literals")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(use_cxx11)

use_cxx11()


# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )



# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------


# DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# Cairo

FIND_PATH(CAIRO_INCLUDE_DIRS cairo.h
    $ENV{CAIRODIR}/include
    $ENV{CAIRODIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/include
    /usr/include
    /usr/include/cairo # MT
    /usr/local/include/cairo #brew, manual
    /sw/include # Fink
    /opt/local/include # DarwinPorts
    /opt/local/include/cairo # DarwinPorts # MT
    /opt/csw/include # Blastwave
    /opt/include/cairo
    /usr/X11R6/include/cairo
)


find_library(CAIRO_LIBRAIRIES
             NAMES cairo
             PATHS /usr/lib
                  /usr/local/lib
                 ENV CAIRO_ROOT
                 ENV LD_LIBRARY_PATH
                 ENV LIBRARY_PATH
             PATH_SUFFIXES cairo
            )

SET(CAIRO_FOUND 0)
IF(CAIRO_LIBRAIRIES AND CAIRO_INCLUDE_DIRS)
  SET(CAIRO_FOUND 1)
ENDIF(CAIRO_LIBRAIRIES AND CAIRO_INCLUDE_DIRS)


IF(CAIRO_FOUND)
    INCLUDE_DIRECTORIES(${CAIRO_INCLUDE_DIRS})
    message(STATUS "Cairo found in INCLUDE=" ${CAIRO_INCLUDE_DIRS})
    message(STATUS "Cairo found in LIBRARY=" ${CAIRO_LIBRAIRIES})
ELSE(CAIRO_FOUND)
    message(FATAL_ERROR "cairo not found. Check the cmake variables associated to this package or disable it." )
ENDIF(CAIRO_FOUND)



SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  ADD_DEFINITIONS(${BOOST_DEFINITIONS} -DBOOST_ALL_NO_LIB)
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(GTVDep ${GTVDep} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY} Boost::program_options)
   SET(GTVLinInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

SET(SRC_COMP
 ${PROJECT_SOURCE_DIR}/src/BasicVectoImageExporter.cpp  
)

SET(SOURCE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/src/BezierCurve.h
   ${PROJECT_SOURCE_DIR}/src/BasicVectoImageExporter.h
   ${PROJECT_SOURCE_DIR}/src/BezierTriangle2.h
   ${PROJECT_SOURCE_DIR}/src/BreadthFirstVisitorWithParent.h
   ${PROJECT_SOURCE_DIR}/src/BreadthFirstVisitorWithParent.ih
   ${PROJECT_SOURCE_DIR}/src/CairoViewer.h
   ${PROJECT_SOURCE_DIR}/src/ImageConnecter.h
   ${PROJECT_SOURCE_DIR}/src/ImageTVRegularization.h
   ${PROJECT_SOURCE_DIR}/src/ImageTriangulation.h
)



add_subdirectory(bin)
add_subdirectory(tests)
