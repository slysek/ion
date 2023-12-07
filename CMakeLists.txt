cmake_minimum_required(VERSION 3.26)
project(ion)

set(CMAKE_CXX_STANDARD 17)

add_executable(ion main.cpp
        export_lammps.cpp
        export_lammps.h)
