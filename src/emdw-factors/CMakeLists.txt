cmake_minimum_required(VERSION 3.13)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE
    affinemvg.hpp
    affinemvg.cc
    categorical.hpp
    mvg.hpp
    conditionaldirichlet.hpp
    dirichletset.hpp
    dirichletset.tcc
    conditionalpolya.hpp
    conditionalpolya.tcc
    dirichlet.hpp
    dirichlet.tcc
    discretetable.hpp
    discretetable.tcc
    factor.cc
    factor.hpp
    factor.tcc
    factoroperator.cc
    factoroperator.hpp
    factoroperator.tcc
    #nonlinmvg.hpp
    #nonlinmvg.cc
    polya.hpp
    polya.tcc
    sqrtmvg.cc
    sqrtmvg.hpp
    treefactor.cc
    treefactor.hpp
)
