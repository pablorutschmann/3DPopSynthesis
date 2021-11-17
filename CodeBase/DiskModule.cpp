//
//  DiskModule.cpp
//  DiskClass
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#include "DiskModule.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <map>

using namespace std;

const bool RadiativeCooling = false;  // this should be ON only if computing cooling timescale, otherwise should be OFF in standard simulations and an input TTemp should be used


DiskModel::DiskModel() {
}

DiskModel::DiskModel(string input_address, string output_address, std::map<std::string, double> Options) {
    /*
    Initialize disk

    INPUTS
    - input folder address
    - output folder address
    - dictionary of options loaded from the option file
    */

    InputAddress = input_address;
    OutputAddress = output_address;

    // Assign options
    MP = Options["MP"];
    RP = Options["RP"];
    RCavity = Options["RCavity"];
    TDisp = Options["TDisp"];
    TTemp = Options["TTemp"];
    TRefilling = Options["TRefilling"];
    GasDrop = Options["GasDrop"];
    DustDrop = Options["DustDrop"];
    Refilling = Options["Refilling"];
    TempDrop = Options["TempDrop"];


    SetDisk();

    // constants with mu = 2.37 and gamma = 1.4 in the code units
    SigmaBoltz = 9.39e-13 / MP;
    G = 3.69e8 * MP / (RP * RP * RP);
    KbMuMp = 709.14 / RP;
    Alpha = 0.004;
    Gamma = 7. / 5.;
    Cv = KbMuMp / (Gamma - 1);

    cout << "SigmaBoltz" << '\t' << SigmaBoltz << '\n';
    cout << "Cv" << '\t' << '\t' << Cv << '\n';
    cout << "G" << '\t' << '\t' << G << '\n';
    cout << "KbMuMp" << '\t' << '\t' << KbMuMp << '\n';
}


void DiskModel::SetDisk() {
    /*-- SETTING DISK PROFILES READING FILES --*/

    int i = 0;
    ifstream InputFile;
    double index;

    InputFile.open(OutputAddress + "/restart/disk.txt");

    /*-- reading from file. Dustbar is the dust profile without satellite accretion, it is used for refilling calculation --*/
    while (InputFile >> index >> R[i] >> Dr[i] >> SigmaGas[i] >> SigmaDust[i] >> SigmaDustBar[i] >> Temp[i] >> Area[i]
                     >> OmegaK[i] >> SigmaExponent[i] >> TempExponent[i] >> Opacity[i] >> WMF[i] >> SWMF[i]) {
        if (RCavity > R[i]) {
            SigmaGas[i] = 0;
            SigmaDust[i] = 0;
            SigmaDustBar[i] = 0;
        }
        Opacity[i] = ComputeOpacity(i);
        i++;
    }

    InputFile.close();

    if (i == 0)          // if no disk is found in the restart folder, then load it from the input folder
    {
        InputFile.open(InputAddress + "/disk.txt");
        while (InputFile >> index >> R[i] >> Dr[i] >> SigmaGas[i] >> SigmaDust[i] >> SigmaDustBar[i] >> Temp[i]
                         >> Area[i] >> OmegaK[i] >> SigmaExponent[i] >> TempExponent[i] >> Opacity[i] >> WMF[i]
                         >> SWMF[i]) {
            if (RCavity > R[i]) {
                SigmaGas[i] = 0;
                SigmaDust[i] = 0;
                SigmaDustBar[i] = 0;
            }
            Opacity[i] = ComputeOpacity(i);
            i++;

        }
    }

    Length = i;
}


double DiskModel::ComputeCs(int i) {
    /*-- SOUND SPEED --*/
    return sqrt(KbMuMp * Temp[i]);
}


double DiskModel::ComputeH(int i) {
    /*-- SCALE HEIGHT --*/
    return ComputeCs(i) / OmegaK[i];
}

double DiskModel::ComputeOpacity(int i) {
    /*-- COMPUTE OPACITY FROM TABLES --*/
    double P, T, k;

    T = Temp[i];
    P = KbMuMp * T * SigmaGas[i] / ComputeH(i) / sqrt(2. * M_PI);       // this is in Jupiter units
    P /= 3.67e-6;        // this is in cgs

    double logP = log10(P), logT = log10(T), logk;

    if (logT < (0.03 * logP + 3.12)) logk = 0.738 * logT - 1.277;
    else if (logT < (0.0281 * logP + 3.19)) logk = -42.98 * logT + 1.312 * logP + 135.1;
    else if (logT < (0.03 * logP + 3.28)) logk = 4.063 * logT - 15.013;
    else if (logT < (0.00832 * logP + 3.41)) logk = -18.48 * logT + 0.676 * logP + 58.93;
    else if (logT < (0.015 * logP + 3.7)) logk = 2.905 * logT + 0.498 * logP - 13.995;
    else if (logT < (0.04 * logP + 3.91)) logk = 10.19 * logT + 0.382 * logP - 40.936;
    else if (logT < (0.28 * logP + 3.69)) logk = -3.36 * logT + 0.928 * logP + 12.026;
    else {
        if ((logk < (3.586 * logT - 16.85)) && (logT < 4)) logk = 3.586 * logT - 16.85;
        else if (logT < 2.9) logk = 0.738 * logT - 1.277;
        else logk = -0.48;
    }

    k = pow(10, logk);    // in cgs unit
    return k * 8481.;      // in Jupiter unit
}

double DiskModel::ComputeTmin(int i) {
    return 5780. * sqrt(RP / R[i]) * 0.7520883995742579 / 100;
    // return 0
}

int DiskModel::ComputeIceLine() {
// Find index of IceLine in disk
    int index = 0;
    if (Temp[0] < TIce) {
        index = 0;
    }
    else {
        for (int i = Length - 1; i > 0; i--) {
            if (Temp[i] < TIce) {
                index = i;
            } else break;
        }
    }
    return index;
}

double DiskModel::ComputeMaxMass(int i) {
    /*-- MAX DUST MASS FOR EACH CELL --*/
    return Area[i] * SigmaDust[i];
}


void DiskModel::DiskEvolution(double dt) {
    /*-- TIME EVOLUTION OF THE DISK --*/

    int i;
    double DiskFactor, TempFactor;

    // exponential factors to be implemented
    DiskFactor = exp(-dt / TDisp);
    TempFactor = exp(-dt / TTemp);

    for (i = 0; i < Length; i++) {
        double sigma_gas = SigmaGas[i] / 2.;
        if (GasDrop) SigmaGas[i] *= DiskFactor;
        if (DustDrop) {
            SigmaDust[i] *= DiskFactor;
            SigmaDustBar[i] *= DiskFactor;
        }
        sigma_gas += SigmaGas[i] / 2.;      // average gas density before and after if RadiativeCooling is ON
        if ((TempDrop) && (Temp[i] > ComputeTmin(i))) {
            if (RadiativeCooling) {
                double tau = sqrt(2. * M_PI) * SigmaGas[i] * ComputeOpacity(i);
                double T = Temp[i];
                double Tmin = ComputeTmin(i);
                double TempGrad =
                        -SigmaBoltz / Cv / sigma_gas * (T * T * T * T - Tmin * Tmin * Tmin * Tmin) / (tau + 1 / tau);
                /*
                cout << "i = " << i << " over " << Length << '\n';
                cout << "dt = " << dt << '\n';
                cout << "T = " << T << '\n';
                cout << "Tmin = " << Tmin << '\n';
                cout << "tau = " << tau << '\n';
                cout << "T gradient = " << TempGrad << "\n\n";
                */


                Temp[i] = T + TempGrad * dt;
                if (Temp[i] < ComputeTmin(i)) { Temp[i] = ComputeTmin(i); }
            } else {
                double temp = ComputeTmin(i) + (Temp[i] - ComputeTmin(i)) * TempFactor;
                Temp[i] = temp;
            }
        }
        if ((GasDrop) || (TempDrop)) Opacity[i] = ComputeOpacity(i);
    }
    IceLineID = ComputeIceLine();
}


void DiskModel::DiskRefilling(double dt) {
    /*-- COMPUTE REFILLING --*/

    if (Refilling) {
        double rate;

        for (int i = 0; i < Length; i++) {
            if (dt > TRefilling) SigmaDust[i] = SigmaDustBar[i];
            else {
                rate = (SigmaDustBar[i] - SigmaDust[i]) / TRefilling;
                SigmaDust[i] += (rate * dt);
            }
        }
    }
}


double DiskModel::GasMass() {
    /*-- COMPUTE TOTAL GAS MASS --*/

    float sum = 0.;
    for (int i = 0; i < Length; i++) {
        sum += (SigmaGas[i] * Area[i]);
    }
    return sum;
}


double DiskModel::DustMass() {
    /*-- COMPUTE TOTAL DUST MASS --*/

    float sum = 0.;
    for (int i = 0; i < Length; i++) {
        sum += (SigmaDust[i] * Area[i]);
    }
    return sum;
}


double DiskModel::DustBarMass() {
    /*-- COMPUTE DUST "BAR" MASS --*/

    float sum = 0.;
    for (int i = 0; i < Length; i++) {
        sum += (SigmaDustBar[i] * Area[i]);
    }
    return sum;
}

