cmake_minimum_required(VERSION 3.13)
project(Levenshtein)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

set(CMAKE_CXX_STANDARD 11)

add_executable(Levenshtein main.cpp MyBarrier.cpp MyBarrier.h)