cmake_minimum_required(VERSION 3.11)

project(emdw)

message(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "Build  directory: ${CMAKE_BINARY_DIR}")
message(STATUS "Configuring ${PROJECT_TITLE} ...")

#define CMAKE_MODULE_PATH to use some of our macros and functions...
#and include some of our functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(MacroOutOfSourceBuild)
include(FetchContentHelpers)
include(ShowCmakeVariables)

macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

# Required boost support
find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

# Required gLinear support
find_package(gLinear REQUIRED)
include_directories(${GLINEAR_INCLUDE_DIRS})

# Required gLinear support
find_package(prlite REQUIRED)
include_directories(${PRLITE_INCLUDE_DIRS})

# Required patrecII support
#find_package(PatrecII REQUIRED)
#include_directories(${PATRECII_INCLUDE_DIRS})

include(CheckLibraryExists)
include(CheckIncludeFiles)

#****************************************************************
# Example: Add docopt to our project
#****************************************************************
#FetchContent_gitdefaults(docopt
#        https://github.com/docopt/docopt.cpp.git 3dd23e3) #3dd23e3 Dec 2018
#add_subdirectory(${docopt_SOURCE_DIR} ${docopt_BINARY_DIR})
#include_directories(${docopt_SOURCE_DIR} ${docopt_BINARY_DIR})

#****************************************************************
#Example: Add pkassert to our project - a bit more hassle
# (1) no CMAkeList.txt, create one to specify our targets
# (2) note the files are one directory deeper: /src
#****************************************************************
#FetchContent_gitdefaults(ppk_assert
#        https://github.com/gpakosz/PPK_ASSERT.git f64625c)
#
#FILE(WRITE ${ppk_assert_SOURCE_DIR}/src/CMakeLists.txt "
#  cmake_minimum_required(VERSION 3.11)
#  set(CMAKE_INCLUDE_CURRENT_DIR ON)
#  add_library(ppk_assert SHARED ppk_assert.cpp)
#  target_link_libraries(ppk_assert)
#")
#
#add_subdirectory(${ppk_assert_SOURCE_DIR}/src ${ppk_assert_BINARY_DIR}/src)
#include_directories(${ppk_assert_SOURCE_DIR}/src ${ppk_assert_BINARY_DIR}/src)


#****************************************************************
# Add GoogleTest to our project as a shared library
#****************************************************************
#FetchContent_gitdefaults(googletest
#	https://github.com/google/googletest.git release-1.8.0)

#add_subdirectory(${googletest_SOURCE_DIR}/googlemock ${googletest_BINARY_DIR}/googlemock)


#****************************************************************
# Add xtl and xtensor to our project (header only projects)
#****************************************************************
#FetchContent_gitdefaults(xtl https://github.com/QuantStack/xtl.git 0.6.5 )
#FetchContent_gitdefaults(xtensor https://github.com/QuantStack/xtensor.git 0.20.8 )

#include_directories(${xtensor_SOURCE_DIR}/include
#	${xtl_SOURCE_DIR}/include)


#****************************************************************
# Our main project targets
#****************************************************************
add_subdirectory(src)
