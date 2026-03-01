cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    lbp_cg.cc
    lbp_cg.hpp
    lbu2_cg.cc
    lbu2_cg.hpp
    lbu_cg.cc
    lbu_cg.hpp
    messagequeue.cc
    messagequeue.hpp
    lbu_sweep.cc
    lbu_sweep.hpp
)
