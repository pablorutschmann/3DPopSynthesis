//
//  main.cpp
//
//  Created by Marco Cilibrasi on 20/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#include "SatelliteModule.hpp"                        //file with the definition of the classes and functions
#include "DiskModule.hpp"
#include "EvolutionModule.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h> 
#include <chrono>
#include <ctime>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
    /*-- READ ARGUMENTS --*/

    string InputAddress = argv[1];
    string OutputAddress = argv[2];
    string HistoryAddress = argv[3];
    ofstream HistoryFile;

    cout << '\n' << '\n';

    /*-- INITIALIZE EVOLUTION --*/

    EvolutionModel evo(InputAddress, OutputAddress);
    
    cout << "Input = " << evo.InputAddress << '\n';
    cout << "Output = " << evo.OutputAddress << '\n' << '\n';
    cout << "NSatellites = " << evo.NSatellites << '\n';


    cout << '\n';

    /*-- UPDATE HISTORY FILE --*/

    HistoryFile.open(HistoryAddress, ios_base::app);

    auto start = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start);
    cout << "started computation at " << ctime(&start_time);

    HistoryFile << InputAddress << " started at " << ctime(&start_time) << '\n';
    HistoryFile.close();

    /*-- START SIMULATION --*/

    evo.Simulation();

    /*-- FINAL OPERATIONS --*/

    HistoryFile.open(HistoryAddress, ios_base::app);

    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);

    double RunTime = elapsed_seconds.count();

    cout << "finished computation at " << ctime(&end_time)
              << "elapsed time: " << RunTime << " s, " << RunTime/60 << " min, " << RunTime/3600 << " hours\n\n";

    HistoryFile << '\t' << InputAddress << " finished at " << ctime(&end_time) << " after " << RunTime << " s, ";
    HistoryFile << RunTime / 60 << " min, " << RunTime / 3600 << " hours\n";
    HistoryFile.close();

    return 0;
}
