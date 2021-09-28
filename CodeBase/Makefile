all: ../Executable/3DPopSyn
3DPopSyn: ../Executable/main.o SatelliteModule.o DiskModule.o EvolutionModule.o NBodyModule.o
	g++ -O2 -std=c++11 main.o SatelliteModule.o DiskModule.o EvolutionModule.o NBodyModule.o -o 3DPopSyn
main.o: main.cpp SatelliteModule.hpp CodeBase/DiskModule.hpp CodeBase/EvolutionModule.hpp
	g++ -O2 -std=c++11 -c main.cpp -o main.o
SatelliteModule.o : SatelliteModule.cpp SatelliteModule.hpp
	g++ -O2 -std=c++11 -c SatelliteModule.cpp -o SatelliteModule.o
DiskModule.o: DiskModule.cpp CodeBase/DiskModule.hpp
	g++ -O2 -std=c++11 -c DiskModule.cpp -o DiskModule.o
EvolutionModule.o: EvolutionModule.cpp CodeBase/EvolutionModule.hpp CodeBase/DiskModule.hpp SatelliteModule.hpp NBodyModule.hpp
	g++ -O2 -std=c++11 -c EvolutionModule.cpp -o EvolutionModule.o
NBodyModule.o:	NBodyModule.cpp NBodyModule.hpp
	g++ -O2 -std=c++11 -c NBodyModule.cpp -o NBodyModule.o
