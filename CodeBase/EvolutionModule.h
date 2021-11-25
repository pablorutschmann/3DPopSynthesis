//
//  EvolutionModule.h
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#ifndef SIMULATION_EVOLUTIONMODULE_H
#define SIMULATION_EVOLUTIONMODULE_H

//#include "DiskModule.h"
//#include "SatelliteModule.h"
//#include "NBodyModule.h"
//#include <cmath>
//#include <fstream>
//#include <string>
//#include <iostream>
//#include <stdio.h>
//#include <map>
//#include <sys/stat.h>
//#include <stdlib.h>
//#include <chrono>
//#include <ctime>
//#include <time.h>

#include "DiskModule.h"
#include "SatelliteModule.h"
#include "NBodyModule.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <map>
#include <sys/stat.h>
#include <chrono>
#include <ctime>
#include <ctime>
#include <cstdlib>

const int TotalNumberSatellites = 1100;
using namespace std;


class EvolutionModel {

public:

    EvolutionModel();

    EvolutionModel(string input_address, string output_address);

    DiskModel Disk;
    SatelliteModel Satellites[TotalNumberSatellites];
    SatelliteModel SatellitesBackUp[TotalNumberSatellites];

    string InputAddress, OutputAddress, MigrationType;
    double Time, MaxTime;
    int NEmbryos;
    int NPlanetesimals;
    int NSatellites;
    double R_min, R_max;
    double EmbryoInitMass, EmbryoRho;
    double Spacing;
    double InitMass, Rho;
    double DtMax, GlobalDt;
    double RDistruction;
    double AccCoeff, FeedRadius, ThresholdMass;
    int PebbleAccretion;
    double StokesNumber, PebbleFlux;
    double DiskPrecision, RotationFraction;
    double UpdateInterval, UpdateTime;
    double SaveInterval;
    double Snapshots;
    double MaxRunTime;
    int SaveIndex;
    int TraceWM;
    int Sublimation;
    double Tsubli;

    double TimeStopFormation;

    double r_prev;
    double RHill_prev;

    int TotalSubticks;
    int CloseSubticks;

    std::map<std::string, double> Options;

    // Functions

    void Simulation();

    void SetOptions();

    void SetParameters();

    void SatelliteInitialization();

    void CreateSatellite(int index, bool type = 0);

    void CreateSnapshot();

    void WriteSnapshot(string FolderName, bool header);

    void ComputeParameters(int index);

    void SortSatellites();

    void CheckCollision(int index);

    bool CheckInvalidity(int index);

    void CheckAndCreate();

    void Accretion(int index, double dt);

    void AccretionPebble(int index, double dt);

    int Tick();

    int SubTick();

    void I(int i, double factor);

    void K(int i, double factor);

    void KEncounter(int id_group, double factor);

    double InitialTimeStep(int id_group);

    double HPC(int id_group, double dt);

    double ComputeK(int i, int j, int der);

    int ComputeIDGroup(int whichcond);

    int ActiveSatellites();

    void DestroySatellite(int index, int code, int index2);

    void PrintFlags(bool all);

    double Energy();


};

string ZeroPadNumber(int num, int length);

double Kij(double y);

double dKijdy(double y);

double d2Kijdy2(double y);

#endif //SIMULATION_EVOLUTIONMODULE_H
