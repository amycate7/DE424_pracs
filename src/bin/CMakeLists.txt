cmake_minimum_required(VERSION 3.13)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(example example.cc)
target_link_libraries(example emdw)

add_executable(week1 week1.cc)
target_link_libraries(week1 emdw)

add_executable(week2_Q1 week2_Q1.cc)
target_link_libraries(week2_Q1 emdw)

add_executable(week2_Q2 week2_Q2.cc)
target_link_libraries(week2_Q2 emdw)

add_executable(week2_Q3 week2_Q3.cc)
target_link_libraries(week2_Q3 emdw)

add_executable(week3_Q1 week3_Q1.cc)
target_link_libraries(week3_Q1 emdw)

add_executable(week3_Q2 week3_Q2.cc)
target_link_libraries(week3_Q2 emdw)

add_executable(week4_Q1 week4_Q1.cc)
target_link_libraries(week4_Q1 emdw)

add_executable(week4_Q2 week4_Q2.cc)
target_link_libraries(week4_Q2 emdw)

add_executable(week5_Q1 week5_Q1.cc)
target_link_libraries(week5_Q1 emdw)

add_executable(week5_Q2 week5_Q2.cc)
target_link_libraries(week5_Q2 emdw)