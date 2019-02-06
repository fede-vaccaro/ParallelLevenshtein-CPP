cmake_minimum_required(VERSION 3.13)
project(Levenshtein)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

set(CMAKE_CXX_STANDARD 11)

add_executable(Levenshtein main.cpp BoostBarrier.cpp BoostBarrier.h Worker.cpp Worker.h ind.cpp ind.h Worker2.cpp Worker2.h  Worker3.cpp Worker3.h  uint.h ConcurrentQueue.cpp ConcurrentQueue.h)