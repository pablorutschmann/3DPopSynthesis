//
//  DiskModule.hpp
//  DiskClass
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#ifndef DiskModule_hpp
#define DiskModule_hpp

#include <stdio.h>
#include <string>
#include <map>

using namespace std;

const int DiskLength = 1000;

// Defining class for disks

class DiskModel
{
public:
    
    DiskModel();
    DiskModel(string input_address, string output_address, std::map<std::string, double> Options);
    
    string InputAddress, OutputAddress;


    double MP, RP;
    double DustToGas;
    double RCavity;
    double TDisp, TTemp, TRefilling;
    int GasDrop, DustDrop, TempDrop, Refilling;
    
    double G, Cv, SigmaBoltz, KbMuMp, Alpha, Gamma;
    
    int Length;
    double R[DiskLength];
    double Dr[DiskLength];
    double SigmaGas[DiskLength];
    double SigmaDust[DiskLength], SigmaDustBackUp[DiskLength];
    double SigmaDustBar[DiskLength];
    double Temp[DiskLength];
    double Area[DiskLength];
    double OmegaK[DiskLength];
    double SigmaExponent[DiskLength];
    double TempExponent[DiskLength];
    double Opacity[DiskLength];
    
    // Setting functions
    void SetDisk();

    // Scalar functions
    double GasMass();
    double DustMass();
    double DustBarMass();
    
    
    // Vector functions
    double ComputeCs(int i);
    double ComputeH(int i);
    double ComputeMaxMass(int i);
    double ComputeOpacity(int i);
    double ComputeTmin(int i);

    // Evolution functions
    void DiskEvolution(double dt);
    void DiskRefilling(double dt);

    
};

#endif /* DiskModule_hpp */
