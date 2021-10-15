//
//  EvolutionModule.hpp
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#ifndef EvolutionModule_hpp
#define EvolutionModule_hpp

#include <stdio.h>
#include <map>
#include "DiskModule.hpp"
#include "SatelliteModule.hpp"

const int TotalNumberSatellites = 30;
using namespace std;


class EvolutionModel
{

public:
    
    EvolutionModel();
    EvolutionModel(string input_address, string output_address);
    
    DiskModel Disk;
    SatelliteModel Satellites[TotalNumberSatellites];
    SatelliteModel SatellitesBackUp[TotalNumberSatellites];
    
    string InputAddress, OutputAddress, MigrationType;
    double Time, MaxTime;
    int NSatellites;
    double R_min, R_max;
    double InitMass, Rho;
    double DtMax, GlobalDt;
    double RDistruction;
    double AccCoeff, FeedRadius, ThresholdMass;
    double DiskPrecision, RotationFraction;
    double UpdateInterval, UpdateTime;
    double SaveInterval;
    double Snapshots;
    double MaxRunTime;
    int SaveIndex;

    double TimeStopFormation;

    int TotalSubticks;
    int CloseSubticks;
    
    std::map<std::string, double> Options;
    
    // Functions

    void Simulation();
    
    void SetOptions();
    void SetParameters();

    void SatelliteInitialization();
    void CreateSatellite(int index);
    void CreateSnapshot();
    void WriteSnapshot(string FolderName, bool header);
    void ComputeParameters(int index);
    void SortSatellites();
    void CheckCollision(int index);
    void CheckAndCreate();
    void Accretion(int index, double dt);

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

#endif /* EvolutionModule_hpp */
