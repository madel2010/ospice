PROJECT(blockklu)

cmake_minimum_required(VERSION 2.6)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

SET(CMAKE_CXX_FLAGS "-O0 -g -fopenmp -Wall")
SET(CMAKE_C_FLAGS "-O0 -g -fopenmp -Wall")
SET(LIBS gomp)

SET(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/Include 
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/AMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/BTF/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CCOLAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/COLAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CSparse/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CXSparse/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/LDL/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/UMFPACK/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/SuiteSparse_config
		  #the SymbolicC++ library
		  ${PROJECT_SOURCE_DIR}/SymbolicC++3-3.35/include)


set(LIBSPARSE_LIBS
	"${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/Lib/libklu.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/UMFPACK/Lib/libumfpack.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/AMD/Lib/libamd.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Lib/libcholmod.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/COLAMD/Lib/libcolamd.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/CCOLAMD/Lib/libccolamd.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/CAMD/Lib/libcamd.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/BTF/Lib/libbtf.a" 
)

#######START: Find Lapack #####################
find_package(PkgConfig)

find_library(HAVE_LAPACK lapack)
if(HAVE_LAPACK)
    SET(LIBS ${LIBS} lapack)
else()
    ##TRY USING PKGCONFIG
    pkg_check_modules(LAPACK lapack)
    if(LAPACK_FOUND)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LAPACK_CFLAGS}")
        set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${LAPACK_CFLAGS}")
        set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
        include_directories(${LAPACK_INCLUDEDIR})
        
        link_directories(${LAPACK_LIBDIR})
    else()
	unset(LAPACK_FOUND CACHE)
	message(FATAL_ERROR "Can not find lapack library:")
    endif()
	
endif()

find_library(HAVE_BLAS blas)
if(HAVE_BLAS)
    SET(LIBS ${LIBS} blas)
else()
    ##try using pckg_config
    pkg_check_modules(BLAS blas)
    if(BLAS_FOUND)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${BLAS_CFLAGS}")
        set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${BLAS_CFLAGS}")
        set(LIBS ${LIBS} ${BLAS_LIBRARIES})
        include_directories(${BLAS_INCLUDEDIR})
        link_directories(${BLAS_LIBDIR})
    else()
	unset(LAPACK_FOUND CACHE)
	message(FATAL_ERROR "Can not find blas library:")
    endif()
    
endif()
######END: Find Lapack #####################


##Find required Boiost libraries#####
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD OFF)
set(Boost_ADDITIONAL_VERSIONS 1.44 1.44.0 1.49 1.49.0 1.50 1.50.0)
FIND_PACKAGE(Boost  COMPONENTS regex)
if(NOT Boost_REGEX_FOUND)
	unset(Boost_REGEX_FOUND CACHE)
	message(FATAL_ERROR "Can not find Boost:regex library")
else()
	set(LIBS ${LIBS} ${Boost_REGEX_LIBRARY})
endif()

######Check boost/lexical_cast.hpp header files exists#############
INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("boost/lexical_cast.hpp" HAVE_BOOST_LEXICAL_HPP)
if(NOT HAVE_BOOST_LEXICAL_HPP)
	unset(HAVE_BOOST_LEXICAL_HPP CACHE)
	message(FATAL_ERROR "boost/lexical_cast.hpp doesnot exist. Please make sure that Boost library is correctly installed")
endif()


CHECK_INCLUDE_FILE_CXX("boost/regex.hpp" HAVE_BOOST_REGEX_HPP)
if(NOT HAVE_BOOST_REGEX_HPP)
	unset(HAVE_BOOST_REGEX_HPP CACHE)
	message(FATAL_ERROR "boost/regex.hpp doesnot exist. Please make sure that Boost library is correctly installed")
endif()

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -O0 -g)

#ADD_SUBDIRECTORY(src)

include_directories(${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/templates)
include_directories(${INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/templates)
link_directories(${PROJECT_SOURCE_DIR}/src)

######START:Make SuiteSParse#############
add_custom_command(
  COMMAND echo "Compiling SuiteSparse Library"
  COMMAND cd SuiteSparse && make -j5 library
  OUTPUT  ${PARENT_SOURCE_DIR}/SuiteSparse/KLU/lib/libklu.a
)

ADD_CUSTOM_TARGET(suitesparse ALL DEPENDS ${PARENT_SOURCE_DIR}/SuiteSparse/KLU/lib/libklu.a)
######END:Make SuiteSParse#############

######START:Make the simulator library#############
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(SymbolicC++3-3.35)
set(LIBS  simulator SymbolicC++ ${LIBS}) 
######END:Make the simulator library#############


#list all source files here
ADD_EXECUTABLE(Ospice main.cpp)
set(LIBS  ${LIBS} ${LIBSPARSE_LIBS}) 
target_link_libraries(Ospice ${LIBS})

#dependecies of the projects
add_dependencies(Ospice suitesparse SymbolicC++ )

#################For DOXYGEN documentaton#######################
  
  FIND_PACKAGE(Doxygen)
  if (DOXYGEN_FOUND)
    message(STATUS "DOXYGEN documentation is enabled")
    #message(FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly")
  
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxygen.in 
                 ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (doc 
  			COMMAND ${DOXYGEN_EXECUTABLE} "${PROJECT_BINARY_DIR}/Doxyfile"
  			SOURCES "${PROJECT_BINARY_DIR}/Doxyfile")
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
  else()
	message(STATUS "DOXYGEN documentation is Not found")
  endif()
