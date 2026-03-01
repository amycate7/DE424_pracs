cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-std=c++17 -O3 ") # -DTESTING_MODE -DDEBUG_MODE -DCHK_IDX -DDEBUG_MODE -pedantic
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ctor-dtor-privacy")

# to enable warnings
add_definitions ("-Wall")

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Add the gLinear include directory
include_directories(${GLINEAR_INCLUDE_DIRS} ${PRLITE_INCLUDE_DIRS} ${PATRECII_INCLUDE_DIRS})

add_library(emdw SHARED "") #populate via add_subdir...()
target_link_libraries(emdw ${GLINEAR_LIBRARIES} ${PRLITE_LIBRARIES} ${BOOST_LIBRARIES} lapack blas ssl crypto z)
# ${PATRECII_LIBRARIES}

add_subdirectory(emdw-base)       #These all assign target_sources to emdw
add_subdirectory(emdw-beliefprop)
add_subdirectory(emdw-clustering)
add_subdirectory(emdw-factors)
add_subdirectory(emdw-graphing)
add_subdirectory(test)
add_subdirectory(utils-dependancyglue)
add_subdirectory(utils-informationtheory)
add_subdirectory(utils-matrix)
add_subdirectory(utils-misc)
add_subdirectory(utils-sigmapoints)
add_subdirectory(utils-vector)
add_subdirectory(vmp)
add_subdirectory(wip)

#Here are some binaries to compile and run
add_subdirectory(bin)
