cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    affinemvg.test.cc
    bronkerbosch.test.cc
    conditionalpolya.test.cc
    dirichlet.test.cc
    dirichletset.test.cc
    discretetable.test.cc
    polya.test.cc
    sigmapoints.test.cc
    sqrtmvg.test.cc
)
