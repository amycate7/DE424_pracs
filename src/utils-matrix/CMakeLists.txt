cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    mathfuncs.cc
    mathfuncs.hpp
    matops.cc
    matops.hpp
    matops.tcc
    v2vtransform.cc
    v2vtransform.hpp
)