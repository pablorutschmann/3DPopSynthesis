//
//  SatelliteModule.h
//  SatelliteClass
//
//  Defining satellite class
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#ifndef SIMULATION_SATELLITEMODULE_H
#define SIMULATION_SATELLITEMODULE_H

#include <stdio.h>
#include <string>

using namespace std;

class SatelliteModel {
public:

    SatelliteModel();

    SatelliteModel(int id, bool type, double mass, double x, double y, double z, double rho, double g, double mp,
                   double rp, double sn, double time);

    int ID, Index;
    bool Type;
    double Mass;
    double WM;
    double SWM;
    double X, Y, Z;
    double Vx, Vy, Vz;
    double Ax, Ay, Az;
    double Adx, Ady, Adz;           // accelerations for Hermite integrator
    double Addx, Addy, Addz;
    double Adddx, Adddy, Adddz;
    double Radius;
    double Rho;
    double E2D_prev;

    double MP;
    double RP;
    double Mu;
    double StokesNumber;
    double Dt;
    double InitTime;
    double FormationTime;
    double Tsubli;
    double kb;

    int N;
    int IClock, KClock;

    bool Active;
    bool AdvanceI;
    int GroupID;

    double A, Ecc, Inc, R, OmegaK, SigmaExp, TempExp, Rfeed, RHill, Cs, H, Eta, SigmaGas, SigmaDust, Twave, Tmig, Tecc, Tinc, Opacity, Temp, P;

    // Functions
    double ComputeRadius();

    double ComputeR();

    double ComputeR2D();

    double ComputeRHill();

    double ComputeV();

    double ComputeTheta();

    void SetInitVel();

    double ComputeEnergy();

    double ComputeEcc();

    double ComputeA();

    double ComputeInc();

    double ComputeP(double alpha);

    double ComputeSublimationRate(double T);

    void SetDt(double global_dt, double rotation_fraction);

    void ComputeTdyn(string migtype, double fraction, double alpha, double gamma, double kbmump, double sigma_boltz);

    double ComputeBmig(string migtype, double alpha, double gamma, double kbmump, double sigma_boltz);

    void ComputeAcc(int MigOption, int EccOption, int IncOption);

    double ComputeE2D();

    void Print(float time, string message);
};

double F(double P);

double G(double P);

double K(double P);

double f(double x, double gamma);

double FJ(double x);

double z(double x, double n);

double GapDepth(double P);


#endif //SIMULATION_SATELLITEMODULE_H
