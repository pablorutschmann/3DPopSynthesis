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
#include <random>

const int TotalNumberSatellites = 501;
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
    int NEmbryos, NPlanetesimals, NSatellites;
    double R_min, R_max;
    double Spacing;
    double EmbryoInitMass, EmbryoRho;
    double InitMass, Rho;
    double DtMax, GlobalDt;
    double RDistruction;
    double AccCoeff, FeedRadius, ThresholdMass;
    double StokesNumber;
    double DiskPrecision, RotationFraction;
    double UpdateInterval, UpdateTime;
    double SaveInterval;
    double Snapshots;
    double MaxRunTime;
    int SaveIndex;
    double SublimationFactor;

    double TimeStopFormation;

    double r_prev;
    double RHill_prev;
    double eff_bar_total;

    int TotalSubticks;
    int CloseSubticks;

    std::map<std::string, double> Options;

    default_random_engine generator;

    // Functions

    void Simulation();

    void SetOptions();

    void SetParameters();

    void SatelliteInitialization();

    void CreateSatellite(int index, bool type = 0);

    double Density_Model(double x);

    double Rejection_Sample();

    void CreateSnapshot();

    void WriteSnapshot(string FolderName, bool header);

    void ComputeParameters(int index);

    void SortSatellites();

    void CheckCollision(int index);

    bool CheckInvalidity(int index);

    void CheckAndCreate();

    void Accretion(int index, double dt);

    double PebbleAccretion(int index, double dt);

    void Sublimation(int index, double dt);

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

double PhysicalRadiusDistribution(double r);


#endif //SIMULATION_EVOLUTIONMODULE_H
