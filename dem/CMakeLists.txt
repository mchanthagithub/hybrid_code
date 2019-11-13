cmake_minimum_required(VERSION 3.10)
project(dem)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")
include_directories(fileIO)
add_executable(dem main.cpp fileIO/outputWriter.cpp fileIO/outputWriter.h Region.cpp Region.h)