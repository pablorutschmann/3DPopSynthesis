#!/bin/zsh

cd CodeBase
make -B
mv 3DPopSyn DiskModule.o EvolutionModule.o main.o NBodyModule.o SatelliteModule.o ../Executable/
cd ..
