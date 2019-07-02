cmake_minimum_required(VERSION 3.13)
project(projet_CPP)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(projet_CPP
        cpp_projet.cpp
        EF2d-base.cpp
        EF2d-base.hpp
        GC.cpp
        GC.hpp
        R2.hpp)
