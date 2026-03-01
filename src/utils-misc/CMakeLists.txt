cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    anytype.hpp
    anytype.tcc
    colourmap.cc
    colourmap.hpp
    combinations.cc
    combinations.hpp
    combinations.tcc
    dtw.hpp
    hashfunctions.hpp
    oddsandsods.cc
    oddsandsods.hpp
)
