cmake_minimum_required(VERSION 3.13)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

target_include_directories(emdw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(emdw PRIVATE)

add_executable(simlda simlda.cc)
target_link_libraries(simlda emdw)

add_executable(emdwlda emdwlda.cc)
target_link_libraries(emdwlda emdw)
