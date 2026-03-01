cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    sortindices.hpp
    sortindices.tcc
    vecset.hpp
    vecset.tcc
)