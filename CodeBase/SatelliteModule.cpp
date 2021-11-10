//
//  SatelliteModule.cpp
//  SatelliteClass
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#include "SatelliteModule.hpp"


#include <iostream>
#include <cmath>
#include <string>

using namespace std;

const int Nmax = pow(2,15); // maximum levels in the adaptive timestep
const bool TypeIIMigration = true;  // take into account gap opening effects


// Defining the functions for satellites

SatelliteModel::SatelliteModel() {
}

SatelliteModel::SatelliteModel(int id, bool type, double mass, double x, double y, double z, double rho, double mu, double time, double sublimationtime)
{
    /*
    Initialize satellite model

    INPUTS
    - satellite ID
    - satellites mass
    - satellite coordinates x,y,z
    - satellite density in unit of Jupiter density
    - mu = G_grav * M_planet
    - time of creation
    */

    ID = id;
    Type = type;
    Mass = mass;
    WM = -1;                    // Water Mass
    SWM = -1;                   // Hydrated Minerals Mass
    X = x;
    Y = y;
    Z = z;
    Rho = rho;
    Mu = mu;                    // G * M
    N = 0;                      // individual timestep index
    Twave = 1e10;               // timescale for migration, eccentricity and inclination damping
    InitTime = time;            // time at which it is initiated
    FormationTime = 0.;         // formation time, updated when the satellites grow to the threshold mass
    Tsubli = sublimationtime;    // Sublimation timescale

    // clockes for individual timesteps, see Saha & Tremaine (1992, 1994)
    IClock = 0;             
    KClock = 0;
    AdvanceI = false;
    
    Active = true;              // whether the satellite is active or destroyed
    
    GroupID = -1;               // group ID for close-encounters
    
    SetInitVel();               // setting initial keplerian velocity
    Radius = ComputeRadius();   // it computes the satellite physical radius

    // acceleration and its derivatives (for Hermite integrator)
    Ax = Ay = Az = 0.;
    Adx = Ady = Adz = 0.;
    Addx = Addy = Addz = 0.;
    Adddx = Adddy = Adddz = 0.;
}

double SatelliteModel::ComputeRadius()
{
    /*
    Compute physical radius
    */
    return pow((Mass / Rho), (1. / 3.));
}


double SatelliteModel::ComputeR2D()
{
    /*
    Compute the 2D distance from the center (R in cylindrical coordinates)
    */
    return sqrt(X*X + Y*Y);
}

double SatelliteModel::ComputeR()
{
    /*
    Compute the 3D distance from the center (r in spherical coordinates)
    */
    return sqrt(X*X + Y*Y + Z*Z);
}

double SatelliteModel::ComputeRHill()
{
    /*
    Compute the Hill radius
    */
    return ComputeA() * pow(Mass / 3, 1./3.);
}


double SatelliteModel::ComputeV()
{
    /*
    Compute absolute value for velocity
    */
    return sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
}


double SatelliteModel::ComputeTheta()
{
    /*
    Compute azimuthal angle (theta in cylindrical coordinates, phi in spherical coordinates)
    */
    return atan2(Y, X);
}


void SatelliteModel::SetInitVel()
{
    /*
    Set the initial keplerian velocity
    */

    double theta = ComputeTheta();
    double r = ComputeR();
    double v = sqrt(Mu / r);
    
    Vx = -sin(theta) * v;
    Vy = cos(theta) * v;
    Vz = 0.;
}


double SatelliteModel::ComputeEnergy()
{
    /*
    Compute total orbital energy
    */
    double r = ComputeR();
    double v2 = Vx*Vx + Vy*Vy + Vz*Vz;
    
    return v2/2. - Mu/r;
}


double SatelliteModel::ComputeEcc()
{
    /*
    Compute eccentricity
    */
    double r = ComputeR();
    double v2 = Vx*Vx + Vy*Vy + Vz*Vz;

    double A = v2 - Mu/r;
    double B = X*Vx + Y*Vy + Z*Vz;
    
    double C2 = (A*X - B*Vx)*(A*X - B*Vx) + (A*Y - B*Vy)*(A*Y - B*Vy) + (A*Z - B*Vz)*(A*Z - B*Vz);
    double C = sqrt(C2);
    
    return C / Mu;
}


double SatelliteModel::ComputeA()
{
    /*
    Compute semi-major axis from orbital energy
    */
    double E = ComputeEnergy();
    return - Mu / 2 / E;
}


double SatelliteModel::ComputeInc()
{
    /*
    Compute inclination
    */
    double h, hx, hy, hz;
    
    hx = Y*Vz - Z*Vy;
    hy = Z*Vx - X*Vz;
    hz = X*Vy - Y*Vx;
    
    h = sqrt(hx*hx + hy*hy + hz*hz);
    
    return acos(hz / h);
}

double SatelliteModel::ComputeP(double alpha)
{
    /*
    Compute the migration parameter for type I/II migration (Crida & Morbidelli 2007)
    */

    return Cs /(OmegaK * A) * (3./4.) * pow(Mass/3, -1./3.) + 50 * alpha / Mass * (Cs / OmegaK / A) * (Cs / OmegaK / A);
}


void SatelliteModel::SetDt(double global_dt, double rotation_fraction)
{
    /*
    Set the timestep of the satellite

    INPUTS
    - global timestep of the simulation (timestep in the inner part of the disk)
    - rotation fraction, i.e. maximum portion of an orbit that the satellite should travel in one time-step. (rotation_fraction = 0.1 means 10 timesteps needed to
                                                            travel a single orbit)
    */

    int N_new;
    double a = ComputeA();
    if(Active == false)
    {
        // maximum dt possible if the satellite is not active
        N = Nmax;
        Dt =  N * global_dt / 2;
    }
    else if(a > 1.)
    {
        double dt = 2 * M_PI * a / sqrt(Mu / a) * rotation_fraction;    // period * rotation_fraction
        if(dt < global_dt)
	    {
	        N = 2;                              // 2 is the minimum index, in order to allow for sub-steps if necessary
	        dt = global_dt;
	    }
        else
        {
            
            double n = log2(dt / global_dt);
            int n1 = (int) n + 1;               // round n to the upper index and 
            N_new = pow(2., n1);
            if(N == 0) N = N_new;
            else if (N_new < N) N = N_new;
            if(N > Nmax) N = Nmax;  
            Dt = N * global_dt / 2;   
        }
    }
    else
    {   
        // minimum dt if the semi-major axis is < 1
        N = 2;
        Dt = global_dt;
    }
    
}

void SatelliteModel::ComputeTdyn(string migtype, double fraction, double alpha, double gamma, double kbmump, double sigma_boltz)
{

    /*
    Compute the dynamical timescales (migration, eccentricity and inclination damping) following Cresswell & Nelson (2008)

    INPUTS
    - migration type identifier
    - rotation_fraction as defined above
    - viscosity alpha
    - adiabatic index gamma
    - kB / (mu * mp) factor
    - Stefan-Boltzmann constant
    */

    double hr, e, i, ihr, ehr, P;
    double b;
    
    hr = H / R;
    i = ComputeInc();
    e = ComputeEcc();

    ehr = e / hr;
    ihr = i / hr;

    P = (1 + pow(ehr / 2.25, 1.2) + pow(ehr / 2.84, 6))/ (1 - pow(ehr / 2.02, 4));

    Tecc = Twave / 0.78 * (1 - 0.14 * ehr * ehr + 0.06 * ehr * ehr * ehr + 0.18 * ehr * ihr * ihr);
    if(abs(Tecc) < (10 * Dt / fraction)) Tecc = 10 * Dt / fraction;
    Tinc = Twave / 0.544 * (1 - 0.3 * ihr * ihr + 0.24 * ihr * ihr * ihr + 0.14 * ehr * ehr * ihr);
    if(abs(Tinc) < (10 * Dt / fraction)) Tinc = 10 * Dt / fraction;
    
    b = ComputeBmig(migtype, alpha, gamma, kbmump, sigma_boltz);

    if(b != 0)
    {
        Tmig = Twave / b / hr / hr * (P + P / abs(P) * (0.07 * ihr + 0.085 * ihr * ihr * ihr * ihr - 0.08 * ehr * ihr * ihr));
        if(abs(Tmig) < (10 * Dt / fraction))
	{
	  if(Tmig == 0) Tmig = 10 * Dt / fraction;
	  else Tmig = 10 * Dt / fraction * Tmig / abs(Tmig);
	  cout << "Satellite " << ID << " has limited migration\n";
	}
	
    }
    else
    {
        Tmig = 1e300;
    }
    
}



double SatelliteModel::ComputeBmig(string migtype, double alpha, double gamma, double kbmump, double sigma_boltz)
{   

    /*
    Compute the bI parameter for type I migration together with the bII parameter
    */

    double bI, bII, b;

    if(migtype == "Tanaka")       // Tanaka et al. (2002)
    {
        bI = 2.7 - 1.1 * SigmaExp;
    }
    else if(migtype == "Dangelo")  // D'Angelo & Lubow (2010)
    {
        bI = 1.36 - 0.62 * SigmaExp - 0.43 * TempExp;
    }
    else if(migtype == "Jimenez_lin")   // linearized version from Jimenez & Masset (2017)
    {
        double q = Mass, h = H/R;
        double Alpha = -SigmaExp, Beta = -TempExp;
        double rho = SigmaGas / H / sqrt(2. * M_PI);
    
        double chi = 16. * (gamma - 1) * sigma_boltz * Temp * Temp * Temp / (3. * rho * rho * Opacity * kbmump);
        double chi_C = R * R * h * h * OmegaK;
    
        // bL
        double bL = -(2.34 - 0.1 * Alpha + 1.5 * Beta) * f(chi / chi_C, gamma);

        bI = bL + (0.46 - 0.96 * Alpha + 1.8 * Beta) / gamma;
        bI = -bI;
    }
    else if(migtype == "Paardekooper")    // Paardekooper et al. (2010); Paardekooper, Baruteau & Kley (2011)
    {
        double b1, b2, b3, b4, b5;
        double csi = 1.1 / pow(gamma, 0.25) * sqrt(Mass * A / H);
        double nu = alpha * Cs * H;
        double Px = sqrt(OmegaK * A*A * csi*csi*csi / (2 * M_PI * nu));
        double Pv = 2. * Px / 3.;
        double p = -SigmaExp, q = -TempExp;

	    b1 = 2.5 + 1.7 * q - 0.1 * p; 
        b2 = - 1.1 * F(Pv) * G(Px) * (1.5 - p);
        b3 = - 0.7 * (1 - K(Pv)) * (1.5 - p);
        b4 = - 7.9 * (q - (gamma - 1)*p) / gamma * F(Pv) * F(Px) * sqrt(G(Pv)*G(Px));
        b5 = - (2.2 - 1.4/gamma) * (q - (gamma - 1)*p) * sqrt((1 - K(Pv)) * (1 - K(Px)));
        bI = (b1 + b2 + b3 + b4 + b5) * 2 / gamma;

    }
    else if(migtype == "Jimenez")      // full version from Jimenez & Masset (2017)
    {
        double q = Mass, h = H/R;
        double Alpha = -SigmaExp, Beta = -TempExp;
        double rho = SigmaGas / H / sqrt(2. * M_PI);
    
        double chi = 16. * (gamma - 1) * sigma_boltz * Temp * Temp * Temp / (3. * rho * rho * Opacity * kbmump);
        double chi_C = R * R * h * h * OmegaK;
    
        // bL
        double bL = -(2.34 - 0.1 * Alpha + 1.5 * Beta) * f(chi / chi_C, gamma);
        
        // bV_CR
        double h1 = h / sqrt(gamma);
        double xs = R * (1.05 * sqrt(q/h1) + 3.4 * pow(q, (7/3)) * pow(h1, -6.)) / (1 + 2 * pow(q, 2.) * pow(h1, -6.));
    
        double nu = alpha * Cs * H;
        double zv = R * nu / (OmegaK * xs * xs * xs);
    
        double eb = 1. / (1 + 30. * h * zv);
    
        double bV_lin = (0.976 - 0.640 * Alpha) / gamma;
    
        double FV = 8. * M_PI / 3. * zv * FJ(zv);
        double bV_UHD = 3./4. * (3./2. - Alpha) * pow((xs / R), 4.) * (h/q) * (h/q);
        double bV_HD = FV * bV_UHD;
    
        double bV = eb * bV_HD + (1. - eb) * bV_lin;

        // bS_CR
    
        double zx = R * chi / (OmegaK * xs * xs * xs);
        double ev = 1. / (1. + (6. * h * zv) * (6. * h * zv));
        double ex = 1. / (1. + 15. * h * zx);
    
        double csi = (Beta - 0.4 * Alpha - 0.64);
    
        double bS_lin = 0.8 * csi / gamma;
    
        double fs = 1.2 * min(1., 1.4 * sqrt(zx)) * min(1., 1.8 * sqrt(zv));
        double bS_UHD = 3.3 * csi * pow((xs / R), 4.) * (h / q) * (h / q);
        double bS_HD = fs * bS_UHD;
    
        double bS = ev * ex * bS_HD + (1. - ev * ex) * bS_lin;

        // bT_CR
    
        double bT_lin = 1.0 * Beta / gamma;
    
        double FT = 1.2 * min(1., 1.8 * sqrt(zv));
        double bT_UHD = 0.73 * Beta * pow((xs / R), 4.) * (h/q) * (h/q);
        double bT_HD = FT * bT_UHD;
    
        double bT = ev * bT_HD + (1. - ev) * bT_lin;

        // bVCT_CR
    
        double bVCT = 4. * M_PI * csi / gamma * pow((xs / R), 4.) * (h/q) * (h/q) * eb * zv * (zv * FJ(zv) - zx * FJ(zx)) / (zv - zx);
    
        // b
    
        double bC = bV + bS + bT + bVCT;
        bI = bL + bC;
        bI = -bI;
    }
    else bI = 0;
    
    // bII

    if(TypeIIMigration)
    {
        bI *= GapDepth(P);

        double B = 4 * M_PI * A * A * SigmaGas;
        double FII = 1 + Mass / B;
            
        double factor = 3 * (2 + SigmaExp + TempExp);
            
        bII = factor / FII * (Cs * Cs * Cs * Cs * alpha) / (OmegaK * OmegaK * OmegaK * OmegaK * A * A * A * A * A * A * SigmaGas * Mass);

        b = z(1/P, 30.) * bI + (1 - z(1/P, 30.)) * bII;
    }
    else b = bI;

    return b;
}


void SatelliteModel::ComputeAcc(int MigOption, int EccOption, int IncOption)
{
    /*
    Compute the acceleration components due to migration, eccentricity and inclination damping (Cresswell & Nelson 2008)
    
    INPUTS
    - migration, eccentricity damping and inclination damping switch (1 if active, 0 if not)
    */

    // Migration
    if(MigOption == 1)
    {
        Ax -= Vx / Tmig;
        Ay -= Vy / Tmig;
        Az -= Vz / Tmig;
    }

    // Eccentricity
    if(EccOption == 1)
    {
        double a, r2, rv;
        rv = X*Vx + Y*Vy + Z*Vz;
        r2 = X*X + Y*Y + Z*Z;
        a = -2 * rv / r2 / Tecc;

        Ax += a * X;
        Ay += a * Y;
        Az += a * Z;
    }

    // Inclination
    if(IncOption == 1)
    {
        Az -= Vz / Tinc;
    }
}

void SatelliteModel::UpdateWM(double dt) {
    /*
    Update the Water Mass according to exponential decay with timescale TSublimation,
    also update the total mass of the satellite

    INPUTS
    - the timestep used for the current position and velocity update
    */

    double factor = exp(- dt / Tsubli);
//    cout << "CHECK" << '\n';
//    cout << Tsubli << '\n';
//    cout << dt << '\n';
//    cout << factor << '\n';
    double WM_diff = WM * (1 - factor);
    Mass -= WM_diff;
    WM *= factor;
}


void SatelliteModel::Print(float time, string message)
{
    cout << "\n\nSatellite " << ID << " at time " << time << '\n';
    cout << message << '\n';
    cout << '\n';
    cout << "ID = " << ID << '\n';
    cout << "Group ID = " << GroupID << '\n';
    cout << "Disk index = " << Index << '\n';
    cout << "Mass = " << Mass << '\n';
    cout << "WM = " << WM << '\n';
    cout << "SWM = " << SWM << '\n';
    cout << "Physical r = " << ComputeR2D() << '\n';
    cout << "Disk r = " << R << '\n';
    cout << "X = " << X << '\n';
    cout << "Y = " << Y << '\n';
    cout << "Z = " << Z << '\n';
    cout << "Vx = " << Vx << '\n';
    cout << "Vy = " << Vy << '\n';
    cout << "Vz = " << Vz << '\n';
    cout << "ecc = " << ComputeEcc() << '\n';
    cout << "inc = " << ComputeInc() << '\n';
    cout << "a = " << ComputeA() << '\n';
    cout << "Period = " << sqrt((2 * M_PI) * (2 * M_PI) * ComputeA() * ComputeA() * ComputeA() / Mu) << '\n';
    cout << "dt = " << Dt << '\n';
    cout << "N = " << N << '\n';
    cout << "Active = " << Active << '\n';
    cout << '\n' << '\n';
}


/*-- SOME FUNCTIONS USED FOR MIGRATION --*/

double F(double P)
{
    return 1 / (1 + (P * P / (1.3 * 1.3)));
}

double G(double P)
{
    if(P < sqrt((8. / 45. / M_PI))) return 16. / 25. * pow((45. * M_PI / 8.), 3./4.) * pow(P, 3./2.);
    else return 1. - 9. / 25. * pow((8. / 45. / M_PI), 4./3.) * pow(P, -8./3.);
}

double K(double P)
{
    if(P < sqrt(28. / 45. / M_PI)) return 16. / 25. * pow((45. * M_PI / 28.), 3./4.) * pow(P, 1.5);
    else return 1 - 9./25. * pow((28. / 45. / M_PI), 4./3.) * pow(P, -8./3.);
}

double f(double x, double gamma)
{
    return (sqrt(x / 2) + 1/gamma) / (sqrt(x/2) + 1);
}


double FJ(double x)
{
    if(x < 4./9.) return 1 - sqrt(x);
    else return 4. / (27. * x);
}

double z(double x, double n)
{
    return 1 / (1 + pow(x, n));
}



double GapDepth(double P)
{   
    if(P>2.4646) return 1 - exp(-pow(P, 0.75) / 3);
    else return (P-0.541)/4;
}