cmake_minimum_required(VERSION 3.12)

project(hesse_curve LANGUAGES C)
add_executable(${PROJECT_NAME} main.c)
target_compile_features(${PROJECT_NAME} PRIVATE c_std_99)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake/)

find_package(gmp REQUIRED)
target_link_libraries(${PROJECT_NAME} gmp)
