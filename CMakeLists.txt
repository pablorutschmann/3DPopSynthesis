cmake_minimum_required(VERSION 3.20)
project(3DPopSynthesis)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)

add_executable(3DPopSynthesis CodeBase/main.cpp CodeBase/EvolutionModule.cpp CodeBase/EvolutionModule.h CodeBase/DiskModule.cpp CodeBase/DiskModule.h CodeBase/SatelliteModule.cpp CodeBase/SatelliteModule.h CodeBase/NBodyModule.cpp CodeBase/NBodyModule.h)
