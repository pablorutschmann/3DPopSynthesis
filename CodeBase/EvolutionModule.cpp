//
//  EvolutionModule.cpp
//  SatelliteClass
//
//  Created by Marco Cilibrasi on 21/02/2019.
//  Copyright © 2019 Marco Cilibrasi. All rights reserved.
//

#include "EvolutionModule.hpp"
#include "DiskModule.hpp"
#include "NBodyModule.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>
#include <time.h>
#include <stdlib.h>

using namespace std;

const double MigrationFactor = 1;
const double FormationPerc = 0.1;
const double GalileanMass = 2.53e-5;        // Europa mass
const double MaxInclination = 0.01;
const int NH = 3; // coefficient to define r_crit
const double eta = 0.02; // first coefficient for HPC timestep
const double etaS = 0.01; // second coefficient for HPC timestep


const bool flags = false;


EvolutionModel::EvolutionModel() {
}

EvolutionModel::EvolutionModel(string input_address, string output_address) {
    /*-- INITIALIZE EVOLUTION --*/

    InputAddress = input_address;
    OutputAddress = output_address;

    TotalSubticks = 0;
    CloseSubticks = 0;

    SetOptions();

    // initialize random seed for satellite formation
    int factor = (int) (Options["DustToGas"] * 1e6 + Options["TRefilling"]);
    srand(time(0) + factor); // We need this to have real random numbers and not pseudo-random numbers.

    Disk = DiskModel(InputAddress, OutputAddress, Options);
    GlobalDt = 2 * M_PI * sqrt(RDistruction * RDistruction * RDistruction / Disk.G / Disk.MP) *
               RotationFraction;   // minimum timestep (inner disk)

    SetParameters();

    SatelliteInitialization();
    SortSatellites();
    UpdateInterval = Disk.TDisp * DiskPrecision;     // set the disk update interval


    CreateSnapshot();

    cout << "Initial time = " << Time << '\n' << '\n';
    cout << "Max number of satellites = " << TotalNumberSatellites << '\n';
    cout << "GlobalDt = " << GlobalDt << '\n';
    cout << "Formation probability = " << FormationPerc << '\n';
    cout << "NH = " << NH << '\n';
    cout << "Eta = " << eta << '\n';
    cout << "EtaS = " << etaS << '\n';
    cout << "Migration factor = " << MigrationFactor << '\n' << '\n';
    cout << "Tracing of WM: " << TraceWM << '\n';


}


void EvolutionModel::Simulation() {
    /*-- START EVOLUTION --*/

    double dt;
    int N_backups = 0;     // count how many backups were necessary for close-encounters

    auto start = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start);

    while (Time < MaxTime) {
        auto end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        time_t end_time = chrono::system_clock::to_time_t(end);

        double RunTime = elapsed_seconds.count();

        if ((RunTime / 60. / 60.) >= MaxRunTime) {
            cout << "\n##################################\n";
            cout << "######   RUN TIME IS OVER   ######\n";
            cout << "##################################\n\n";
            break;
        }


        // set global time-step as the maximum time-step
        int N = ActiveSatellites();
        if (N == 0) dt = DtMax;
        else dt = Satellites[NSatellites - 1].Dt;
        if (flags) cout << "dt = " << dt << '\n';


        /*-- BACK-UP SAVING IN CASE OF TICK RESTART --*/

        for (int i = 0; i < NSatellites; i++) {
            SatellitesBackUp[i] = Satellites[i];
        }
        for (int i = 0; i < Disk.Length; i++) {
            Disk.SigmaDustBackUp[i] = Disk.SigmaDust[i];
        }


        /*-- RUN TICKS, see Saha & Tremaine (1992, 1994) for Tick and SubTick formalism --*/

        if (N != 0) {
            int flag;
            while (true) {
                flag = Tick();
                if (flag == 1)    // if back-up needed for close-encounters restart from back-up
                {
                    cout << "\n ---->  Back-up needed at time " << Time << '\n';
                    N_backups += 1;
                    for (int i = 0; i < NSatellites; i++) {
                        Satellites[i] = SatellitesBackUp[i];
                    }
                    for (int i = 0; i < Disk.Length; i++) {
                        Disk.SigmaDust[i] = Disk.SigmaDustBackUp[i];
                    }
                    SortSatellites();
                } else break;
            }
        }


        /*-- FORMATION OF NEW SATELLITES --*/

        if (Time <= TimeStopFormation) CheckAndCreate();

        /*-- DUST REFILLING --*/
        Disk.DiskRefilling(dt);

        Time += dt;

        /*-- UPDATE DISK QUANTITIES --*/
        if ((Time - UpdateTime) >= UpdateInterval) {
            Disk.DiskEvolution(Time - UpdateTime);
            UpdateTime = Time;
        }

        /*-- UPDATE SATELLITE PARAMETERS --*/
        for (int i = 0; i < NSatellites; i++) {
            ComputeParameters(i);
        }
        SortSatellites();

        /*-- CREATE SNAPSHOT --*/
        if ((Time - SaveIndex * SaveInterval) >= SaveInterval) {
            SaveIndex++;
            CreateSnapshot();
        }
    }
    //SaveIndex ++;
    WriteSnapshot(OutputAddress + "/restart", false);

    cout << "Total Subticks = " << TotalSubticks << '\n';
    cout << "Close encounter = " << CloseSubticks << '\n';
    cout << "Ratio = " << (double) CloseSubticks / TotalSubticks << '\n';
    cout << "N of back-ups = " << N_backups << '\n';

    cout << "\nFinal time = " << Time << '\n' << '\n';
}


void EvolutionModel::SetOptions() {
    /*-- SET INITIAL SIMULATION OPTIONS --*/

    ifstream InputFile;
    InputFile.open(InputAddress + "/options.txt");

    string key;
    double value;

    cout << '\n' << '\n';
    while (InputFile >> key >> value) {
        Options[key] = value;
        cout << key << '\t' << Options[key] << '\n';
    }

    InputFile.close();
    NEmbryos = Options["NEmbryos"];
    NPlanetesimals = Options["NPlanetesimals"];
    NSatellites = NEmbryos + NPlanetesimals;
    R_min = Options["R_min"];
    R_max = Options["R_max"];
    EmbryoInitMass = Options["EmbryoInitMass"];
    EmbryoRho = Options["EmbryoRho"];
    Spacing = Options["Spacing"];
    InitMass = Options["InitMass"];
    Rho = Options["Rho"];
    FeedRadius = Options["FeedRadius"];
    ThresholdMass = Options["ThresholdMass"];
    RDistruction = Options["RDistruction"];
    MaxTime = Options["MaxTime"];
    DtMax = Options["DtMax"];
    Snapshots = Options["Snapshots"];
    SaveInterval = Options["SaveInterval"];
    AccCoeff = Options["AccCoeff"];
    RotationFraction = Options["RotationFraction"];
    DiskPrecision = Options["DiskPrecision"];
    MaxRunTime = Options["MaxRunTime"];
    TraceWM = Options["TraceWM"];
    Sublimation = Options["Sublimation"];
    Tsubli = Options["Tsubli"];
    if (Options.find("MigrationType") != Options.end()) {
        int MigIndex = Options["MigrationType"];
        if (MigIndex == 0) MigrationType = "Tanaka";
        else if (MigIndex == 1) MigrationType = "Dangelo";
        else if (MigIndex == 2) MigrationType = "Paardekooper";
        else if (MigIndex == 3) MigrationType = "Jimenez_lin";
        else if (MigIndex == 4) MigrationType = "Jimenez";
        else MigrationType = "Dangelo";
    } else MigrationType = "Jimenez";


    cout << "Migration type\t" << MigrationType << '\n';

}


void EvolutionModel::SetParameters() {
    /*-- SET SIMULATION PARAMETERS --*/

    std::map<std::string, double> Parameters;
    ifstream InputFile;
    InputFile.open(OutputAddress + "/restart/parameters.txt");

    string key;
    double value;
    int i = 0;

    while (InputFile >> key >> value)    // look for parameters in the restart folders
    {
        Parameters[key] = value;
        i++;
    }

    InputFile.close();

    if (i == 0)    // if NO parameters are found in the restart folder create the satellite files
    {

        ofstream OutputFile;

        OutputFile.open(OutputAddress + "/satellite_list.txt");
        OutputFile << "#ID\tType\tinit_time\tmass\twm\tswm\tr2D\ttheta\tx\ty\tz\tinit_temp\n";
        OutputFile.close();

        OutputFile.open(OutputAddress + "/lost_satellites.txt");
        OutputFile
                << "#ID\ttime\tType\tmass\twm\tswm\tr\tx\ty\tz\tvx\tvy\tvz\ta\tecc\tinc\tformation_time\tcollision\tcollision_index\n";
        OutputFile.close();

        OutputFile.open(OutputAddress + "/collisions.txt");
        OutputFile
                << "#time\tID1\tType1\tmass1\twm1\tswm1\tx1\ty1\tz1\tvx1\tvy1\tvz1\tID2\tType2\tmass2\twm2\tswm2\tx2\ty2\tz2\tvx2\tvy2\tvz2\tID\tType\tmass\twm\tswm\tx\ty\tz\tvx\tvy\tvz\n";
        OutputFile.close();

        Time = 0.;
        UpdateTime = 0.;
        SaveIndex = 0;
        Disk.IceLineID = Disk.ComputeIceLine();

        double DustMass = Disk.DustBarMass();
        TimeStopFormation = Time + Disk.TDisp * log(DustMass * FormationPerc / InitMass);
        cout << "\nTimeStopFormation = " << TimeStopFormation << '\n';
    } else {
        Time = Parameters["Time"];
        UpdateTime = Parameters["UpdateTime"];
        SaveIndex = Parameters["SaveIndex"];
        TimeStopFormation = Parameters["TimeStopFormation"];
        cout << "\nTimeStopFormation = " << TimeStopFormation << '\n';
        Disk.IceLineID = Parameters["IceLineID"];
    }
}


void EvolutionModel::SatelliteInitialization() {
    /*-- INITIALIZE SATELLITES FROM RESTART FILES OR FROM SCRATCH --*/

    int ID, N;
    double mass, rho, wm, swm, x, y, z, vx, vy, vz, init_time, form_time, a, e, inc, dt, p;
    bool type;
    ifstream InputFile;
    cout << OutputAddress + "/restart/satellites.txt'\n";
    InputFile.open(OutputAddress + "/restart/satellites.txt");

    cout << "Searching satellites in " << OutputAddress << "/restart/satellites.txt\n";

    int i = 0;

    while (InputFile >> ID >> type >> mass >> wm >> swm >> x >> y >> z >> vx >> vy >> vz >> a >> e >> inc >> N >> dt >> init_time >> form_time >> p)    // search satellites in the restart file
    {
        if (type == 1) {
            rho = EmbryoRho;
        } else {
            rho = Rho;
        }
        Satellites[i] = SatelliteModel(ID, type, mass, x, y, z, rho, Disk.G * Disk.MP, init_time, Tsubli);
        Satellites[i].Vx = vx;
        Satellites[i].Vy = vy;
        Satellites[i].Vz = vz;
        Satellites[i].FormationTime = form_time;
        Satellites[i].N = N;
        Satellites[i].Dt = N * GlobalDt / 2;
        Satellites[i].WM = wm;
        Satellites[i].SWM = swm;
        ComputeParameters(i);

        i++;
    }
    InputFile.close();

    cout << i << " satellites detected\n";

    if (i == 0)   // if satellites are not found, create them
    {
        string folder_name = OutputAddress + "/restart";
        mkdir(folder_name.c_str(), ACCESSPERMS);

        int j = 0;
        for (j = 0; j < NSatellites; j++) {
            if (j < NEmbryos) {
                CreateSatellite(j, 1);
            } else {
                CreateSatellite(j, 0);
            }

            ComputeParameters(j);
        }
    }

    for (i = 0; i < NSatellites; i++) {
        SatellitesBackUp[i] = Satellites[i];
    }

}


void EvolutionModel::CreateSatellite(int index, bool type) {
    /*-- CREATE ONE SATELLITE IN SPECIFIC POSITION --*/

    double r, expr, theta, z, mass, rho;
    int i;
    int ID = 0;

    for (i = 0; i < (NSatellites); i++)   // find next ID
    {
        if (Satellites[i].ID > ID) {
            ID = Satellites[i].ID;
        }
    }
    ID++;

    theta = (rand() % 10000) / 10000. * 2 * M_PI;


    if (ID == 1 and type == 1) {
        r = R_min + 0.1 * R_min;
        z = (2 * (rand() % 10000) / 10000. - 1) * MaxInclination * r;
        r_prev = r;
        Satellites[index] = SatelliteModel(ID, type, EmbryoInitMass, r * cos(theta), r * sin(theta), z, EmbryoRho,
                                           Disk.G * Disk.MP, Time, Tsubli);
        ComputeParameters(index);
        RHill_prev = Satellites[index].ComputeRHill();

    }
    if (ID > 1 and type == 1) {
        r = r_prev + Spacing * RHill_prev +
            ((-0.1 * Spacing * RHill_prev) + (rand() % 10) / 10. * (0.2 * Spacing * RHill_prev));
        z = (2 * (rand() % 10000) / 10000. - 1) * MaxInclination * r;
        r_prev = r;
        Satellites[index] = SatelliteModel(ID, type, EmbryoInitMass, r * cos(theta), r * sin(theta), z, EmbryoRho,
                                           Disk.G * Disk.MP, Time, Tsubli);
        ComputeParameters(index);
        RHill_prev = Satellites[index].ComputeRHill();
    }


    if (type == 0) {
        bool invalid = true;
        while (invalid == true) {
            expr = log10(R_min) + (rand() % 10000) / 10000. * (log10(R_max) - log10(R_min));
            r = pow(10, expr);
            z = (2 * (rand() % 10000) / 10000. - 1) * MaxInclination * r;
            Satellites[index] = SatelliteModel(ID, type, InitMass, r * cos(theta), r * sin(theta), z, Rho,
                                               Disk.G * Disk.MP, Time, Tsubli);
            ComputeParameters(index);
            invalid = CheckInvalidity(index);
        }
    }

    ofstream OutputFile;
    OutputFile.open(OutputAddress + "/satellite_list.txt", ios_base::app);
    OutputFile << Satellites[index].ID << '\t' << Satellites[index].Type << '\t' << Satellites[index].InitTime << '\t'
               << Satellites[index].Mass << '\t' << Satellites[index].WM << '\t' << Satellites[index].SWM << '\t' << Satellites[index].ComputeR2D()
               << '\t' << Satellites[index].ComputeTheta() << '\t' << Satellites[index].X << '\t' << Satellites[index].Y
               << '\t' << Satellites[index].Z << '\t' << Disk.Temp[Satellites[index].Index] << '\n';
    OutputFile.close();

    Satellites[index].Print(Time, "");

    SatellitesBackUp[index] = Satellites[index];
}


void EvolutionModel::CreateSnapshot() {
    /*-- CREATE A SNAPSHOT --*/

    int index = SaveIndex;
    string FolderName, StringIndex;

    StringIndex = ZeroPadNumber(index, 4);
    FolderName = OutputAddress + "/Snapshot_" + StringIndex;
    mkdir(FolderName.c_str(), ACCESSPERMS);

    if (Snapshots) {
        WriteSnapshot(FolderName, true);
    }
    WriteSnapshot(OutputAddress + "/restart", false);

    cout << "\nSnapshot " << index << " created at time " << Time << "\n\n";
}


void EvolutionModel::WriteSnapshot(string FolderName, bool header) {
    /*-- WRITE SNAPSHOT WITH OR WITHOUT HEADER, no header is for restart files --*/

    ofstream OutputFile;

    OutputFile.open(FolderName + "/disk.txt");

    if (header)
        OutputFile
                << "#index\tr\tdr\tSigmaGas\tSigmaDust\tSigmaDustBar\tTemp\tArea\tOmegaK\tSigmaExponent\tTExponent\tOpacity\tWMF\tSWMF\n";

    for (int i = 0; i < Disk.Length; i++) {
        OutputFile << i << '\t' << Disk.R[i] << '\t' << Disk.Dr[i] << '\t' << Disk.SigmaGas[i] << '\t'
                   << Disk.SigmaDust[i] << '\t' << Disk.SigmaDustBar[i] << '\t' << Disk.Temp[i] << '\t' << Disk.Area[i]
                   << '\t' << Disk.OmegaK[i] << '\t' << Disk.SigmaExponent[i] << '\t' << Disk.TempExponent[i] << '\t'
                   << Disk.Opacity[i] << '\t' << Disk.WMF[i] << '\n';
    }

    OutputFile.close();


    OutputFile.open(FolderName + "/satellites.txt");

    if (header)
        OutputFile << "#ID\tType\tM\tWM\tSWM\tx\ty\tz\txv\tvy\tvz\ta\te\ti\tN\tdt\tinit_time\tformation_time\tP\n";

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].Active) {
            OutputFile << Satellites[i].ID << '\t' << Satellites[i].Type << '\t' << Satellites[i].Mass << '\t'
                       << Satellites[i].WM << '\t' << Satellites[i].SWM << '\t' << Satellites[i].X << '\t'
                       << Satellites[i].Y << '\t' << Satellites[i].Z << '\t' << Satellites[i].Vx << '\t'
                       << Satellites[i].Vy << '\t' << Satellites[i].Vz << '\t' << Satellites[i].ComputeA() << '\t'
                       << Satellites[i].ComputeEcc() << '\t' << Satellites[i].ComputeInc() << '\t' << Satellites[i].N
                       << '\t' << Satellites[i].Dt << '\t' << Satellites[i].InitTime << '\t'
                       << Satellites[i].FormationTime << '\t' << Satellites[i].P << '\n';
        }
    }

    OutputFile.close();


    OutputFile.open(FolderName + "/parameters.txt");

    OutputFile << "Time" << '\t' << Time << '\n';
    OutputFile << "UpdateTime" << '\t' << UpdateTime << '\n';
    OutputFile << "SaveIndex" << '\t' << SaveIndex << '\n';
    OutputFile << "TimeStopFormation" << '\t' << TimeStopFormation << '\n';
    OutputFile << "GlobalDt" << '\t' << GlobalDt << '\n';
    OutputFile << "IceLineID" << '\t' << Disk.IceLineID << '\n';
    OutputFile << "IceLineRadius" << '\t' << Disk.R[Disk.IceLineID] << '\n';

    OutputFile.close();


}


void EvolutionModel::ComputeParameters(int index) {
    /*-- UPDATE ALL SATELLITE PARAMETERS --*/

    double position;

    if (Satellites[index].Active) {
        // find satellite index in the disk grid 
        position = abs(Satellites[index].ComputeA() * cos(Satellites[index].ComputeInc()));
        if (position < Disk.R[0]) {
            Satellites[index].Index = 0;
        } else if (position >= Disk.R[Disk.Length - 1]) {
            Satellites[index].Index = Disk.Length - 1;
            //
        } else {
            for (int i = 0; i < Disk.Length - 1; i++) {
                if (position < Disk.R[i + 1]) {
                    Satellites[index].Index = i;
                    break;
                }
            }
        }

        Satellites[index].Radius = Satellites[index].ComputeRadius();
        Satellites[index].A = Satellites[index].ComputeA();
        Satellites[index].Ecc = Satellites[index].ComputeEcc();
        Satellites[index].Inc = Satellites[index].ComputeInc();
        Satellites[index].R = Disk.R[Satellites[index].Index];
        Satellites[index].OmegaK = Disk.OmegaK[Satellites[index].Index];
        Satellites[index].SigmaExp = Disk.SigmaExponent[Satellites[index].Index];
        Satellites[index].TempExp = Disk.TempExponent[Satellites[index].Index];
        Satellites[index].RHill = Satellites[index].ComputeRHill();
        Satellites[index].Rfeed = Satellites[index].RHill * FeedRadius;
        Satellites[index].Cs = Disk.ComputeCs(Satellites[index].Index);
        Satellites[index].H = Disk.ComputeH(Satellites[index].Index);
        Satellites[index].SigmaGas = Disk.SigmaGas[Satellites[index].Index];
        Satellites[index].SetDt(GlobalDt, RotationFraction);
        Satellites[index].Opacity = Disk.Opacity[Satellites[index].Index];
        Satellites[index].Temp = Disk.Temp[Satellites[index].Index];
        Satellites[index].P = Satellites[index].ComputeP(Disk.Alpha);

        if (Satellites[index].WM < 0.) {

            Satellites[index].WM = Disk.WMF[Satellites[index].Index] * Satellites[index].Mass;
        }
        Satellites[index].R = Disk.R[Satellites[index].Index];

        if (Satellites[index].SWM < 0.) {

            Satellites[index].SWM = Disk.SWMF[Satellites[index].Index] * Satellites[index].Mass;
        }
        Satellites[index].R = Disk.R[Satellites[index].Index];


        Satellites[index].SigmaDust = Disk.SigmaDust[Satellites[index].Index];
        if (flags)
            cout << "Satellite " << index << " -> Mass = " << Satellites[index].Mass << ", SigmaGas = "
                 << Satellites[index].SigmaGas << '\n';
        if ((Satellites[index].Mass != 0) && (Satellites[index].SigmaGas != 0)) {
            Satellites[index].Twave = 1 / (Satellites[index].Mass * Satellites[index].SigmaGas * position * position *
                                           Satellites[index].OmegaK) * pow((Satellites[index].H / position), 4);
        } else Satellites[index].Twave = 1e300;

        Satellites[index].Twave /= MigrationFactor;

        if ((Satellites[index].Mass >= GalileanMass) & (Satellites[index].FormationTime == 0)) {
            Satellites[index].FormationTime = Time - Satellites[index].InitTime;
        }
    } else {
        Satellites[index].SetDt(GlobalDt, RotationFraction);
    }

}


void EvolutionModel::SortSatellites() {
    /*-- SORT SATELLITES AND BACK-UPS BY TIMESTEPS --*/

    bool change = true;
    SatelliteModel temp;

    while (change)   // iterate until no change happens
    {
        change = false;
        for (int i = 0; i < NSatellites - 1; i++) {
            if (Satellites[i].N > Satellites[i + 1].N) {
                change = true;

                temp = Satellites[i];
                Satellites[i] = Satellites[i + 1];
                Satellites[i + 1] = temp;

                temp = SatellitesBackUp[i];
                SatellitesBackUp[i] = SatellitesBackUp[i + 1];
                SatellitesBackUp[i + 1] = temp;
            }
        }
    }
}


double EvolutionModel::ComputeK(int i, int j, int der) {
    /*-- COMPUTE K COEFFICIENT AND ITS DERIVATIVES FOR COUPLES OF SATELLITES --*/
    if (i == j) {
        if (der == 0) return 1;
        else return 0;
    }


    double distance, distance2;
    double Xij, Yij, Zij;
    double rcrit;
    rcrit = NH * pow((Satellites[i].Mass + Satellites[j].Mass) / 3, 1. / 3.) * (Satellites[i].A + Satellites[j].A) / 2;

    Xij = Satellites[i].X - Satellites[j].X;
    Yij = Satellites[i].Y - Satellites[j].Y;
    Zij = Satellites[i].Z - Satellites[j].Z;

    distance2 = Xij * Xij + Yij * Yij + Zij * Zij;
    distance = sqrt(distance2);

    double y;
    if (rcrit <= 0) {
        y = 2;
        rcrit = 1.;
    } else y = (distance - 0.1 * rcrit) / (0.9 * rcrit);

    if (der == 0) return Kij(y);
    else if (der == 1) return dKijdy(y) / (0.9 * rcrit);
    else if (der == 2) return d2Kijdy2(y) / ((0.9 * rcrit) * (0.9 * rcrit));
    else return 0;
}


int EvolutionModel::Tick() {
    /*-- COMPUTE A TICK, FLAG = 1 -> RESTART TICK FROM BACK-UP  --*/
    int FinalN = Satellites[NSatellites - 1].N;
    int flag;

    for (int i = 0; i < NSatellites; i++) {
        Satellites[i].KClock = 0.;
        Satellites[i].IClock = 0.;
    }
    while (true) {
        flag = SubTick();
        if (flag == 1) return 1;
        else {
            bool exit = true;
            for (int i = 0; i < NSatellites; i++) {
                if (Satellites[i].KClock != FinalN) {
                    exit = false;
                    break;
                }
            }
            if (exit) break;
        }
    }

    if (flags) { cout << "\nTick is done at time " << Time << "\n\n\n\n"; }

    return 0;
}


int EvolutionModel::ComputeIDGroup(int whichcond) {
    /*-- DETECT GROUPS OF CLOSE ENCOUNTERS BEFORE APPLYING Ks (whichcond says which of the 2 Ks) --*/

    bool is_close_encounter = false;
    double K;

    for (int i = 0; i < NSatellites; i++) {
        bool cond1, cond2, cond;
        cond1 = (Satellites[i].IClock == Satellites[i].KClock);
        cond2 = ((i == 0) || ((Satellites[i].IClock == Satellites[i - 1].IClock) && (Satellites[i - 1].AdvanceI)));
        if (whichcond == 1) cond = cond1;
        else if (whichcond == 2) cond = cond2;

        if (cond) {
            Satellites[i].GroupID = -1;
        }
    }

    if (Options["NBody"] == false) return 0;


    for (int i = 0; i < NSatellites; i++) {
        bool cond1, cond2, cond;
        cond1 = (Satellites[i].IClock == Satellites[i].KClock);
        cond2 = ((i == 0) || ((Satellites[i].IClock == Satellites[i - 1].IClock) && (Satellites[i - 1].AdvanceI)));
        if (whichcond == 1) cond = cond1;
        else if (whichcond == 2) cond = cond2;

        if (cond && Satellites[i].Active) {
            if (Satellites[i].Mass < ThresholdMass) {
                continue;
            } else {

                for (int j = i + 1; j < NSatellites; j++) {
                    if (j == i) continue;
                    else if (Satellites[j].Active == false) continue;
                    else if (Satellites[j].Mass < ThresholdMass) continue;
                    else {
                        K = ComputeK(i, j, 0);
                    }

                    if (K < 1) {
                        is_close_encounter = true;
                        if (Satellites[i].N != Satellites[j].N) {

                            cout << "Different timestep for satellites " << i << " (ID " << Satellites[i].ID << ") and "
                                 << j << " (ID " << Satellites[j].ID << ")" << '\n';

                            double Dt_to_assign = min(Satellites[i].Dt, Satellites[j].Dt);
                            int N_to_assign = min(Satellites[i].N, Satellites[j].N);

                            cout << Satellites[i].N << '\t' << SatellitesBackUp[i].N << '\t' << Satellites[j].N << '\t'
                                 << SatellitesBackUp[j].N << '\n';

                            Satellites[j].Dt = Dt_to_assign;
                            Satellites[j].N = N_to_assign;
                            SatellitesBackUp[j].Dt = Dt_to_assign;
                            SatellitesBackUp[j].N = N_to_assign;

                            Satellites[i].Dt = Dt_to_assign;
                            Satellites[i].N = N_to_assign;
                            SatellitesBackUp[i].Dt = Dt_to_assign;
                            SatellitesBackUp[i].N = N_to_assign;

                            cout << Satellites[i].N << '\t' << SatellitesBackUp[i].N << '\t' << Satellites[j].N << '\t'
                                 << SatellitesBackUp[j].N << '\n';

                            if (Satellites[j].KClock != Satellites[j].IClock) {
                                cout << "Detected back-up necessity for " << i << " and " << j << '\n';

                                return 1;
                            }
                        }
                        if (Satellites[i].GroupID == -1) Satellites[i].GroupID = Satellites[i].ID;
                        Satellites[j].GroupID = Satellites[i].GroupID;
                    }
                }
            }
        }
    }
    SortSatellites();

    if (is_close_encounter) CloseSubticks += 1;
    TotalSubticks += 1;

    return 0;
}


int EvolutionModel::SubTick() {
    /*-- PROCESSING A SUB-TICK --*/

    int SincCheck = 0;
    SincCheck = ComputeIDGroup(1);
    if (SincCheck == 1) return 1;

    // Flags for debugging 
    PrintFlags(true);

    if (flags) {
        for (int i = 0; i < NSatellites; i++)
            cout << Satellites[i].ID << ": K = " << Satellites[i].KClock << '\t' << "I = " << Satellites[i].IClock
                 << '\n';
        cout << '\n';
    }


    // Proceed with the 1st Ks

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].KClock == Satellites[i].IClock) {
            if (Satellites[i].GroupID == -1) K(i, 0.5);
            else if (Satellites[i].GroupID == Satellites[i].ID) KEncounter(Satellites[i].GroupID, 0.5);

            Satellites[i].KClock += Satellites[i].N / 2;
            Satellites[i].AdvanceI = true;
            if (Options["Accretion"]) {
                if (Satellites[i].Active) Accretion(i, Satellites[i].Dt * 0.5);
            }
        }

        if (flags)
            cout << Satellites[i].ID << ": K = " << Satellites[i].KClock << '\t' << "I = " << Satellites[i].IClock
                 << '\n';
    }

    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].Active) && (Satellites[i].AdvanceI)) CheckCollision(i);
    }

    if (flags) cout << '\n';

    // Proceed with the Is

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].AdvanceI) {
            I(i, 1.);

            Satellites[i].IClock += Satellites[i].N;
            Satellites[i].AdvanceI = false;
        }

        if (flags)
            cout << Satellites[i].ID << ": K = " << Satellites[i].KClock << '\t' << "I = " << Satellites[i].IClock
                 << '\n';
    }

    if (flags) cout << '\n';

    // Proceed with the 2nd Ks

    SincCheck = 0;
    SincCheck = ComputeIDGroup(2);
    if (SincCheck == 1) return 1;

    PrintFlags(false);

    for (int i = 0; i < NSatellites; i++) {
        if ((i == 0) || ((Satellites[i].IClock == Satellites[i - 1].IClock) && (Satellites[i - 1].AdvanceI))) {
            if (Satellites[i].GroupID == -1) K(i, 0.5);
            else if (Satellites[i].GroupID == Satellites[i].ID) KEncounter(Satellites[i].GroupID, 0.5);
            Satellites[i].KClock += Satellites[i].N / 2;
            Satellites[i].AdvanceI = true;
            if (Options["Accretion"]) {
                if (Satellites[i].Active) Accretion(i, Satellites[i].Dt * 0.5);
            }
        }
        if (flags)
            cout << Satellites[i].ID << ": K = " << Satellites[i].KClock << '\t' << "I = " << Satellites[i].IClock
                 << '\n';
    }

    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].Active) && (Satellites[i].AdvanceI)) CheckCollision(i);
    }

    if (flags) cout << '\n';

    if (flags) double en = Energy();

    return 0;
}


void EvolutionModel::K(int i, double factor) {
    /*-- COMPUTE KEPLERIAN EVOLUTION --*/

    double dt = Satellites[i].Dt * factor;
    double mu = Disk.G * Disk.MP;
    double *x, *y, *z, *vx, *vy, *vz;
    int check;
    ofstream OutputFile;

    if (Options["Sublimation"] == true) {

        double position = abs(Satellites[i].ComputeA() * cos(Satellites[i].ComputeInc()));
        if (position < Disk.R[Disk.IceLineID]) {
            Satellites[i].UpdateWM(dt);
        }
    }

    x = &(Satellites[i].X);
    y = &(Satellites[i].Y);
    z = &(Satellites[i].Z);
    vx = &(Satellites[i].Vx);
    vy = &(Satellites[i].Vy);
    vz = &(Satellites[i].Vz);

    check = fgfull(x, y, z, vx, vy, vz, dt, mu);

    /* Destroy Satellites with hyperbolic trajectory */
    if (check && Satellites[i].Active) DestroySatellite(i, 2, -1);
}


void EvolutionModel::KEncounter(int id_group, double factor) {
    /*-- COMPUTE KEPLERIAN + ENCOUNTER EVOLUTION --*/

    double dt_total = 0.;

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].GroupID == id_group) {
            dt_total = Satellites[i].Dt * factor;
            break;
        }
    }
    double dt = InitialTimeStep(id_group);

    while (true) {
        if (dt_total <= 0) break;
        if (dt >= dt_total) dt = dt_total;
        if (flags) cout << "HPC dt = " << dt << " with dt_total = " << dt_total << '\n';
        dt_total -= dt;
        dt = HPC(id_group, dt);
        for (int i = 0; i < NSatellites; i++) {
            if (Satellites[i].GroupID == id_group) {
                if (Satellites[i].Active) CheckCollision(i);
            }
        }
    }
}


double EvolutionModel::InitialTimeStep(int id_group) {
    /*-- COMPUTE INITIAL TIME-STEP FOR HERMITE INTEGRATION OF A GIVEN GROUP --*/

    double dt = 0;

    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].GroupID == id_group) && (Satellites[i].Active)) {
            double ax, ay, az, adx, ady, adz;
            double r2, Ri, Xi, Yi, Zi, Vxi, Vyi, Vzi, Ri3, Ri5, RiVi;
            double mu = Disk.G * Disk.MP;

            Xi = Satellites[i].X;
            Yi = Satellites[i].Y;
            Zi = Satellites[i].Z;

            r2 = (Xi * Xi) + (Yi * Yi) + (Zi * Zi);
            Ri = sqrt(r2);
            Ri3 = Ri * Ri * Ri;
            Ri5 = Ri3 * Ri * Ri;

            Vxi = Satellites[i].Vx;
            Vyi = Satellites[i].Vy;
            Vzi = Satellites[i].Vz;

            RiVi = (Xi * Vxi) + (Yi * Vyi) + (Zi * Vzi);

            ax = -mu * Xi / Ri3;
            ay = -mu * Yi / Ri3;
            az = -mu * Zi / Ri3;
            adx = -mu * (Vxi / Ri3 - 3 * Xi * RiVi / Ri5);
            ady = -mu * (Vyi / Ri3 - 3 * Yi * RiVi / Ri5);
            adz = -mu * (Vzi / Ri3 - 3 * Zi * RiVi / Ri5);

            for (int j = 0; j < NSatellites; j++) {
                if ((Satellites[j].GroupID == id_group) && (i != j) && (Satellites[j].Active)) {
                    double r2, Rij, Xij, Yij, Zij, Vxij, Vyij, Vzij, Rij3, Rij5, RijVij;


                    Xij = Satellites[i].X - Satellites[j].X;
                    Yij = Satellites[i].Y - Satellites[j].Y;
                    Zij = Satellites[i].Z - Satellites[j].Z;

                    r2 = (Xij * Xij) + (Yij * Yij) + (Zij * Zij);
                    Rij = sqrt(r2);
                    Rij3 = Rij * Rij * Rij;
                    Rij5 = Rij3 * Rij * Rij;

                    Vxij = Satellites[i].Vx - Satellites[j].Vx;
                    Vyij = Satellites[i].Vy - Satellites[j].Vy;
                    Vzij = Satellites[i].Vz - Satellites[j].Vz;

                    RijVij = (Xij * Vxij) + (Yij * Vyij) + (Zij * Vzij);

                    double factor = 1 - ComputeK(i, j, 0) + Rij * ComputeK(i, j, 1);
                    double factor2 = RijVij * ComputeK(i, j, 2) / (Rij * Rij * Rij);
                    double F = -Disk.G * Satellites[j].Mass;

                    ax += (F * Xij / Rij3) * factor;
                    ay += (F * Yij / Rij3) * factor;
                    az += (F * Zij / Rij3) * factor;


                    adx += (F * ((Vxij / Rij3 - 3 * Xij * RijVij / Rij5) * factor + Xij * factor2));
                    ady += (F * ((Vyij / Rij3 - 3 * Yij * RijVij / Rij5) * factor + Yij * factor2));
                    adz += (F * ((Vzij / Rij3 - 3 * Zij * RijVij / Rij5) * factor + Zij * factor2));

                }
            }


            double dti;
            double a, ad;

            a = sqrt(ax * ax + ay * ay + az * az);
            ad = sqrt(adx * adx + ady * ady + adz * adz);

            dti = etaS * a / ad;

            if (dt == 0) dt = dti;
            else if (dti < dt) dt = dti;
        }
    }

    return dt;
}


double EvolutionModel::HPC(int id_group, double dt) {
    /*-- COMPUTE HERMITE INTEGRATION FOR A GROUP OF PARTICLES AND NEW TIME-STEP --*/




    // If you modify the first for loop you have to modify the second loop accordingly

    double next_dt = 0;

    // calculate acceleration and derivative of acceleration
    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].GroupID == id_group) && Satellites[i].Active) {
            double ax, ay, az, adx, ady, adz;
            double r2, Ri, Xi, Yi, Zi, Vxi, Vyi, Vzi, Ri3, Ri5, RiVi;
            double mu = Disk.G * Disk.MP;


            Xi = Satellites[i].X;
            Yi = Satellites[i].Y;
            Zi = Satellites[i].Z;

            r2 = (Xi * Xi) + (Yi * Yi) + (Zi * Zi);
            Ri = sqrt(r2);
            Ri3 = Ri * Ri * Ri;
            Ri5 = Ri3 * Ri * Ri;

            Vxi = Satellites[i].Vx;
            Vyi = Satellites[i].Vy;
            Vzi = Satellites[i].Vz;

            RiVi = (Xi * Vxi) + (Yi * Vyi) + (Zi * Vzi);

            ax = -mu * Xi / Ri3;
            ay = -mu * Yi / Ri3;
            az = -mu * Zi / Ri3;
            adx = -mu * (Vxi / Ri3 - 3 * Xi * RiVi / Ri5);
            ady = -mu * (Vyi / Ri3 - 3 * Yi * RiVi / Ri5);
            adz = -mu * (Vzi / Ri3 - 3 * Zi * RiVi / Ri5);

            for (int j = 0; j < NSatellites; j++) {
                if ((Satellites[j].GroupID == id_group) && (i != j) && (Satellites[j].Active)) {
                    double r2, Rij, Xij, Yij, Zij, Vxij, Vyij, Vzij, Rij3, Rij5, RijVij;


                    Xij = Satellites[i].X - Satellites[j].X;
                    Yij = Satellites[i].Y - Satellites[j].Y;
                    Zij = Satellites[i].Z - Satellites[j].Z;

                    r2 = (Xij * Xij) + (Yij * Yij) + (Zij * Zij);
                    Rij = sqrt(r2);
                    Rij3 = Rij * Rij * Rij;
                    Rij5 = Rij3 * Rij * Rij;

                    Vxij = Satellites[i].Vx - Satellites[j].Vx;
                    Vyij = Satellites[i].Vy - Satellites[j].Vy;
                    Vzij = Satellites[i].Vz - Satellites[j].Vz;

                    RijVij = (Xij * Vxij) + (Yij * Vyij) + (Zij * Vzij);

                    double factor = 1 - ComputeK(i, j, 0) + Rij * ComputeK(i, j, 1);
                    double factor2 = RijVij * ComputeK(i, j, 2) / (Rij * Rij * Rij);
                    double F = -Disk.G * Satellites[j].Mass;

                    ax += (F * Xij / Rij3) * factor;
                    ay += (F * Yij / Rij3) * factor;
                    az += (F * Zij / Rij3) * factor;


                    adx += (F * ((Vxij / Rij3 - 3 * Xij * RijVij / Rij5) * factor + Xij * factor2));
                    ady += (F * ((Vyij / Rij3 - 3 * Yij * RijVij / Rij5) * factor + Yij * factor2));
                    adz += (F * ((Vzij / Rij3 - 3 * Zij * RijVij / Rij5) * factor + Zij * factor2));

                }
            }

            Satellites[i].Ax = ax;
            Satellites[i].Ay = ay;
            Satellites[i].Az = az;

            Satellites[i].Adx = adx;
            Satellites[i].Ady = ady;
            Satellites[i].Adz = adz;
        }
    }

    // update velocities and position accordingly
    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].GroupID == id_group) && Satellites[i].Active) {
            double ax, ay, az, adx, ady, adz;

            ax = Satellites[i].Ax;
            ay = Satellites[i].Ay;
            az = Satellites[i].Az;

            adx = Satellites[i].Adx;
            ady = Satellites[i].Ady;
            adz = Satellites[i].Adz;

            Satellites[i].X += Satellites[i].Vx * dt + ax * dt * dt / 2 + adx * dt * dt * dt / 6;
            Satellites[i].Y += Satellites[i].Vy * dt + ay * dt * dt / 2 + ady * dt * dt * dt / 6;
            Satellites[i].Z += Satellites[i].Vz * dt + az * dt * dt / 2 + adz * dt * dt * dt / 6;

            Satellites[i].Vx += ax * dt + adx * dt * dt / 2;
            Satellites[i].Vy += ay * dt + ady * dt * dt / 2;
            Satellites[i].Vz += az * dt + adz * dt * dt / 2;
        }
    }

    // interpolate and derive corrections for positions and velocities
    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].GroupID == id_group) && Satellites[i].Active) {
            double ax, ay, az, adx, ady, adz;
            double r2, Ri, Xi, Yi, Zi, Vxi, Vyi, Vzi, Ri3, Ri5, RiVi;
            double mu = Disk.G * Disk.MP;
            double a0x, a0y, a0z, ad0x, ad0y, ad0z;
            double addx, addy, addz, adddx, adddy, adddz;

            a0x = Satellites[i].Ax;
            a0y = Satellites[i].Ay;
            a0z = Satellites[i].Az;

            ad0x = Satellites[i].Adx;
            ad0y = Satellites[i].Ady;
            ad0z = Satellites[i].Adz;

            Xi = Satellites[i].X;
            Yi = Satellites[i].Y;
            Zi = Satellites[i].Z;

            r2 = (Xi * Xi) + (Yi * Yi) + (Zi * Zi);
            Ri = sqrt(r2);
            Ri3 = Ri * Ri * Ri;
            Ri5 = Ri3 * Ri * Ri;

            Vxi = Satellites[i].Vx;
            Vyi = Satellites[i].Vy;
            Vzi = Satellites[i].Vz;

            RiVi = (Xi * Vxi) + (Yi * Vyi) + (Zi * Vzi);

            ax = -mu * Xi / Ri3;
            ay = -mu * Yi / Ri3;
            az = -mu * Zi / Ri3;
            adx = -mu * (Vxi / Ri3 - 3 * Xi * RiVi / Ri5);
            ady = -mu * (Vyi / Ri3 - 3 * Yi * RiVi / Ri5);
            adz = -mu * (Vzi / Ri3 - 3 * Zi * RiVi / Ri5);

            for (int j = 0; j < NSatellites; j++) {
                if ((Satellites[j].GroupID == id_group) && (i != j) && Satellites[j].Active) {
                    double r2, Rij, Xij, Yij, Zij, Vxij, Vyij, Vzij, Rij3, Rij5, RijVij;


                    Xij = Satellites[i].X - Satellites[j].X;
                    Yij = Satellites[i].Y - Satellites[j].Y;
                    Zij = Satellites[i].Z - Satellites[j].Z;

                    r2 = (Xij * Xij) + (Yij * Yij) + (Zij * Zij);
                    Rij = sqrt(r2);
                    Rij3 = Rij * Rij * Rij;
                    Rij5 = Rij3 * Rij * Rij;

                    Vxij = Satellites[i].Vx - Satellites[j].Vx;
                    Vyij = Satellites[i].Vy - Satellites[j].Vy;
                    Vzij = Satellites[i].Vz - Satellites[j].Vz;

                    RijVij = (Xij * Vxij) + (Yij * Vyij) + (Zij * Vzij);

                    double factor = 1 - ComputeK(i, j, 0) + Rij * ComputeK(i, j, 1);
                    double factor2 = RijVij * ComputeK(i, j, 2) / (Rij * Rij * Rij);
                    double F = -Disk.G * Satellites[j].Mass;

                    ax += (F * Xij / Rij3) * factor;
                    ay += (F * Yij / Rij3) * factor;
                    az += (F * Zij / Rij3) * factor;


                    adx += (F * ((Vxij / Rij3 - 3 * Xij * RijVij / Rij5) * factor + Xij * factor2));
                    ady += (F * ((Vyij / Rij3 - 3 * Yij * RijVij / Rij5) * factor + Yij * factor2));
                    adz += (F * ((Vzij / Rij3 - 3 * Zij * RijVij / Rij5) * factor + Zij * factor2));

                }
            }

            addx = (-6 * (a0x - ax) - dt * (4 * ad0x + 2 * adx)) / (dt * dt);
            addy = (-6 * (a0y - ay) - dt * (4 * ad0y + 2 * ady)) / (dt * dt);
            addz = (-6 * (a0z - az) - dt * (4 * ad0z + 2 * adz)) / (dt * dt);

            adddx = (12 * (a0x - ax) + 6 * dt * (ad0x + adx)) / (dt * dt * dt);
            adddy = (12 * (a0y - ay) + 6 * dt * (ad0y + ady)) / (dt * dt * dt);
            adddz = (12 * (a0z - az) + 6 * dt * (ad0z + adz)) / (dt * dt * dt);

            Satellites[i].Addx = addx;
            Satellites[i].Addy = addy;
            Satellites[i].Addz = addz;

            Satellites[i].Adddx = adddx;
            Satellites[i].Adddy = adddy;
            Satellites[i].Adddz = adddz;

            // compute next dt
            double dti;
            double a, ad, add, addd;

            a = sqrt(ax * ax + ay * ay + az * az);
            ad = sqrt(adx * adx + ady * ady + adz * adz);
            add = sqrt(addx * addx + addy * addy + addy * addy);
            addd = sqrt(adddx * adddx + adddy * adddy + adddz * adddz);

            dti = sqrt(eta * (a * add + ad * ad) / (ad * addd + add * add));

            if (next_dt == 0) next_dt = dti;
            else if (dti < next_dt) next_dt = dti;
        }
    }


    //Implement corrections
    for (int i = 0; i < NSatellites; i++) {
        if ((Satellites[i].GroupID == id_group) && Satellites[i].Active) {
            double addx, addy, addz, adddx, adddy, adddz;

            if (Options["Sublimation"] == true) {
                double position = abs(Satellites[i].ComputeA() * cos(Satellites[i].ComputeInc()));
                if (position < Disk.R[Disk.IceLineID]) {
                    Satellites[i].UpdateWM(dt);
                }
            }

            addx = Satellites[i].Addx;
            addy = Satellites[i].Addy;
            addz = Satellites[i].Addz;

            adddx = Satellites[i].Adddx;
            adddy = Satellites[i].Adddy;
            adddz = Satellites[i].Adddz;

            Satellites[i].X += (addx * dt * dt * dt * dt / 24 + adddx * dt * dt * dt * dt * dt / 120);
            Satellites[i].Y += (addy * dt * dt * dt * dt / 24 + adddy * dt * dt * dt * dt * dt / 120);
            Satellites[i].Z += (addz * dt * dt * dt * dt / 24 + adddz * dt * dt * dt * dt * dt / 120);

            Satellites[i].Vx += (addx * dt * dt * dt / 6 + adddx * dt * dt * dt * dt / 24);
            Satellites[i].Vy += (addy * dt * dt * dt / 6 + adddy * dt * dt * dt * dt / 24);
            Satellites[i].Vz += (addz * dt * dt * dt / 6 + adddz * dt * dt * dt * dt / 24);
        }
    }

    return next_dt;
}


void EvolutionModel::I(int i, double factor) {
    /*-- COMPUTE INTERACTION PART OF THE HAMILTONIAN --*/

    double dt = Satellites[i].Dt * factor;

    if (Options["Sublimation"] == true) {
        double position = abs(Satellites[i].ComputeA() * cos(Satellites[i].ComputeInc()));
        if (position < Disk.R[Disk.IceLineID]) {
            Satellites[i].UpdateWM(dt);
        }
    }

    Satellites[i].Ax = 0.;
    Satellites[i].Ay = 0.;
    Satellites[i].Az = 0.;

    Satellites[i].ComputeTdyn(MigrationType, RotationFraction, Disk.Alpha, Disk.Gamma, Disk.KbMuMp, Disk.SigmaBoltz);
    Satellites[i].ComputeAcc(Options["Migration"], Options["EccDamping"], Options["IncDamping"]);

    if (flags)
        cout << "\nSatellite " << i << " -> Tmig = " << Satellites[i].Tmig << ", Twave = " << Satellites[i].Twave
             << '\n';

    Satellites[i].Vx += Satellites[i].Ax * dt;
    Satellites[i].Vy += Satellites[i].Ay * dt;
    Satellites[i].Vz += Satellites[i].Az * dt;


    if (Options["NBody"]) {

        for (int j = i + 1; j < NSatellites; j++) {
            if (Satellites[j].Active == false) continue;
            double Xij, Yij, Zij, r, r2, Force, factor;

            Xij = Satellites[i].X - Satellites[j].X;
            Yij = Satellites[i].Y - Satellites[j].Y;
            Zij = Satellites[i].Z - Satellites[j].Z;
            r2 = Xij * Xij + Yij * Yij + Zij * Zij;
            r = sqrt(r2);

            factor = ComputeK(i, j, 0) - r * ComputeK(i, j, 1);
            Force = Disk.G / (r2 * r) * factor;
            Satellites[i].Vx -= Xij * Force * Satellites[j].Mass * dt;
            Satellites[i].Vy -= Yij * Force * Satellites[j].Mass * dt;
            Satellites[i].Vz -= Zij * Force * Satellites[j].Mass * dt;

            Satellites[j].Vx += Xij * Force * Satellites[i].Mass * dt;
            Satellites[j].Vy += Yij * Force * Satellites[i].Mass * dt;
            Satellites[j].Vz += Zij * Force * Satellites[i].Mass * dt;

        }
    }
}

bool EvolutionModel::CheckInvalidity(int index) {
    /*-- CHECK IF NEWLY CREATED SATELLITE IS INSIDE ANOTHER SATELLITE --*/
    bool flag = false;
    for (int j = 0; j < NSatellites; j++) {
        if (j == index) continue;
        if (Satellites[j].Active == false) continue;

        double R1, R2, r2, r;
        R1 = Satellites[index].RHill;

        r2 = (Satellites[index].X - Satellites[j].X) * (Satellites[index].X - Satellites[j].X) +
             (Satellites[index].Y - Satellites[j].Y) * (Satellites[index].Y - Satellites[j].Y) +
             (Satellites[index].Z - Satellites[j].Z) * (Satellites[index].Z - Satellites[j].Z);
        r = sqrt(r2);

        if (r < R1) {
            flag = true;
            break;
        }
    }
    return flag;
}

void EvolutionModel::CheckCollision(int index) {
    /*-- CHECK FOR COLLISIONS OF SATELLITE index WITH SATELLITES OR WITH CENTRAL PLANET OR EJECTION --*/

    double r, r2;
    ofstream OutputFile;

    r = Satellites[index].ComputeR2D();

    if (Options["NBody"] == false) {
        r = 0.;
    } else if (r < RDistruction) {
        DestroySatellite(index, 0, -1);
    } else {
        for (int j = 0; j < NSatellites; j++) {
            if (j == index) continue;
            if (Satellites[j].Active == false) continue;

            double R1, R2;

            if ((Satellites[index].Mass >= ThresholdMass) && (Satellites[j].Mass >= ThresholdMass)) {
                R1 = Satellites[index].Radius;
                R2 = Satellites[j].Radius;
            } else {
                R1 = Satellites[index].RHill;
                R2 = Satellites[j].RHill;
            }
            r2 = (Satellites[index].X - Satellites[j].X) * (Satellites[index].X - Satellites[j].X) +
                 (Satellites[index].Y - Satellites[j].Y) * (Satellites[index].Y - Satellites[j].Y) +
                 (Satellites[index].Z - Satellites[j].Z) * (Satellites[index].Z - Satellites[j].Z);
            r = sqrt(r2);

            if (r < (R1 + R2)) {
                int save_index, dest_index;
                if (Satellites[index].Mass >= Satellites[j].Mass) {
                    save_index = index;
                    dest_index = j;
                } else {
                    save_index = j;
                    dest_index = index;
                }

                OutputFile.open(OutputAddress + "/collisions.txt", ios_base::app);
                OutputFile << Time << '\t';
                OutputFile << Satellites[save_index].ID << '\t' << Satellites[save_index].Type << '\t'
                           << Satellites[save_index].Mass << '\t' << Satellites[save_index].WM << '\t'
                           << Satellites[save_index].SWM << '\t' << Satellites[save_index].X << '\t'
                           << Satellites[save_index].Y << '\t' << Satellites[save_index].Z << '\t'
                           << Satellites[save_index].Vx << '\t' << Satellites[save_index].Vy << '\t'
                           << Satellites[save_index].Vz << '\t';
                OutputFile << Satellites[dest_index].ID << '\t' << Satellites[dest_index].Type << '\t'
                           << Satellites[dest_index].Mass << '\t' << Satellites[dest_index].WM << '\t'
                           << Satellites[dest_index].SWM << '\t' << Satellites[dest_index].X << '\t'
                           << Satellites[dest_index].Y << '\t' << Satellites[dest_index].Z << '\t'
                           << Satellites[dest_index].Vx << '\t' << Satellites[dest_index].Vy << '\t'
                           << Satellites[dest_index].Vz << '\t';

                double TotalMass = Satellites[save_index].Mass + Satellites[dest_index].Mass;
                Satellites[save_index].X = (Satellites[save_index].Mass * Satellites[save_index].X +
                                            Satellites[dest_index].Mass * Satellites[dest_index].X) / TotalMass;
                Satellites[save_index].Y = (Satellites[save_index].Mass * Satellites[save_index].Y +
                                            Satellites[dest_index].Mass * Satellites[dest_index].Y) / TotalMass;
                Satellites[save_index].Z = (Satellites[save_index].Mass * Satellites[save_index].Z +
                                            Satellites[dest_index].Mass * Satellites[dest_index].Z) / TotalMass;
                Satellites[save_index].Vx = (Satellites[save_index].Mass * Satellites[save_index].Vx +
                                             Satellites[dest_index].Mass * Satellites[dest_index].Vx) / TotalMass;
                Satellites[save_index].Vy = (Satellites[save_index].Mass * Satellites[save_index].Vy +
                                             Satellites[dest_index].Mass * Satellites[dest_index].Vy) / TotalMass;
                Satellites[save_index].Vz = (Satellites[save_index].Mass * Satellites[save_index].Vz +
                                             Satellites[dest_index].Mass * Satellites[dest_index].Vz) / TotalMass;

                if (TraceWM) {
                    double TotalWM = Satellites[save_index].WM + Satellites[dest_index].WM;
                    Satellites[save_index].WM = TotalWM;

                    double TotalSWM = Satellites[save_index].SWM + Satellites[dest_index].SWM;
                    Satellites[save_index].SWM = TotalSWM;
                }

                Satellites[save_index].Mass = TotalMass;
                if (Satellites[save_index].Type == 1 or Satellites[dest_index].Type == 1) {
                    Satellites[save_index].Rho = EmbryoRho;
                    Satellites[save_index].Type = 1;
                }

                OutputFile << Satellites[save_index].ID << '\t' << Satellites[save_index].Type << '\t'
                           << Satellites[save_index].Mass << '\t' << Satellites[save_index].WM << '\t'
                           << Satellites[save_index].SWM << '\t' << Satellites[save_index].X << '\t'
                           << Satellites[save_index].Y << '\t' << Satellites[save_index].Z << '\t'
                           << Satellites[save_index].Vx << '\t' << Satellites[save_index].Vy << '\t'
                           << Satellites[save_index].Vz << '\n';
                OutputFile.close();

                DestroySatellite(dest_index, 1, save_index);

                if (dest_index == index) break;

            }
        }
    }
}


void EvolutionModel::CheckAndCreate() {
    /*-- IF SATELLITES ARE GONE, CREATE NEW ONES --*/

    double dust;

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].Active == false) {
            dust = Disk.DustMass();
            if (FormationPerc * dust >= InitMass) {
                CreateSatellite(i);
                for (int i = 0; i < Disk.Length; i++) {
                    Disk.SigmaDust[i] *= (1 - InitMass / dust);
                }
            } else break;
        }
    }
}


void EvolutionModel::Accretion(int index, double dt) {
    /*-- COMPUTE DUST ACCRETION ONTO SATELLITES AND DUST DEPLETION FROM THE DISK --*/

    double DMmax = 0., Mdot, DM;
    double R = Satellites[index].R;

    for (int i = 0; i < Disk.Length; i++) {
        if (abs(R - Disk.R[i]) <= Satellites[index].Rfeed) {
            DMmax += (Disk.SigmaDust[i] * Disk.Area[i]);
        }
    }

    if (DMmax != 0) {
        Mdot = 2 * sqrt(Satellites[index].Radius / R) * Satellites[index].SigmaDust * R * R *
               sqrt(Satellites[index].Mass) * Satellites[index].OmegaK * AccCoeff;
        DM = Mdot * dt;
        if (DM > DMmax) DM = DMmax;


        Satellites[index].Mass += DM;

        for (int i = 0; i < Disk.Length; i++) {
            if (abs(R - Disk.R[i]) <= Satellites[index].Rfeed) {
                Disk.SigmaDust[i] *= (1 - DM / DMmax);
            }
        }
    }

}


int EvolutionModel::ActiveSatellites() {
    int sum = 0;

    for (int i = 0; i < NSatellites; i++) {
        if (Satellites[i].Active) sum++;
    }

    return sum;

}


void EvolutionModel::DestroySatellite(int index, int code, int index2) {
    /*-- DESTROY SATELLITES AND WRITE ON FILES --*/

    ofstream OutputFile;
    OutputFile.open(OutputAddress + "/lost_satellites.txt", ios_base::app);
    OutputFile << Satellites[index].ID << '\t' << Time << '\t' << Satellites[index].Type << '\t'
               << Satellites[index].Mass << '\t' << Satellites[index].WM << '\t' << Satellites[index].SWM << '\t'
               << Satellites[index].ComputeR2D() << '\t' << Satellites[index].X << '\t' << Satellites[index].Y << '\t'
               << Satellites[index].Z << '\t' << Satellites[index].Vx << '\t' << Satellites[index].Vy << '\t'
               << Satellites[index].Vz << '\t' << Satellites[index].ComputeA() << '\t' << Satellites[index].ComputeEcc()
               << '\t' << Satellites[index].ComputeInc() << '\t' << Satellites[index].FormationTime << '\t' << code
               << '\t' << Satellites[index2].ID << '\n';
    OutputFile.close();
    cout << "Satellite " << Satellites[index].ID << " destroyed because of reason " << code << " with satellite "
         << Satellites[index2].ID << '\n';
    Satellites[index].Active = false;
    Satellites[index].Mass = 0.;
}


void EvolutionModel::PrintFlags(bool all) {
    if (flags) {
        if (all) {
            cout << "\n\n";
            cout << "Time = " << Time << '\n';

            cout << "IDs: ";
            for (int i = 0; i < NSatellites; i++) {
                cout << Satellites[i].ID << '\t';
            }
            cout << '\n';

            cout << "Ns: ";
            for (int i = 0; i < NSatellites; i++) {
                cout << Satellites[i].N << '\t';
            }
            cout << '\n';

            cout << "Dts: ";
            for (int i = 0; i < NSatellites; i++) {
                cout << Satellites[i].Dt << '\t';
            }
            cout << '\n';
        }

        cout << "Group IDs: ";
        for (int i = 0; i < NSatellites; i++) {
            cout << Satellites[i].GroupID << '\t';
        }
        cout << '\n';

        cout << "K matrix: \n";
        for (int i = 0; i < NSatellites; i++) {
            for (int j = 0; j < NSatellites; j++) {
                cout << ComputeK(i, j, 0) << '\t';
            }
            cout << '\n';
        }
        cout << '\n';
    }
}


double EvolutionModel::Energy() {
    double kinetic = 0, potential = 0, energy = 0;
    double mi, mj, G = Disk.G, Mp = Disk.MP;
    double vx, vy, vz, v2;
    double xi, yi, zi, ri;
    double xij, yij, zij, rij;
    for (int i = 0; i < NSatellites; i++) {
        xi = Satellites[i].X;
        yi = Satellites[i].Y;
        zi = Satellites[i].Z;
        ri = sqrt(xi * xi + yi * yi + zi * zi);

        vx = Satellites[i].Vx;
        vy = Satellites[i].Vy;
        vz = Satellites[i].Vz;
        v2 = vx * vx + vy * vy + vz * vz;

        mi = Satellites[i].Mass;

        kinetic += mi * v2 / 2;
        potential -= G * Mp * mi / ri;

        for (int j = i + 1; j < NSatellites; j++) {
            xij = xi - Satellites[j].X;
            yij = yi - Satellites[j].Y;
            zij = zi - Satellites[j].Z;

            mj = Satellites[j].Mass;

            rij = sqrt(xij * xij + yij * yij + zij * zij);

            potential -= G * mi * mj / rij;
        }

        energy = kinetic + potential;
    }

    cout << "Time = " << Time << " -> Energy = " << energy << ", Kinetic = " << kinetic << ", Potential = " << potential
         << "\n";

    return energy;
}


string ZeroPadNumber(int num, int length) {
    string ret;

    // the number is converted to string with the help of stringstream
    ret = to_string(num);

    // Append zero chars
    int str_length = ret.length();
    for (int i = 0; i < length - str_length; i++)
        ret = "0" + ret;
    return ret;
}


double Kij(double y) {
    double den = 3 * y * y - 3 * y + 1;

    if (y <= 0) return 0.;
    else if (y >= 1) return 1.;
    else return (y * y * y) / den;
}


double dKijdy(double y) {
    double den = 3 * y * y - 3 * y + 1;

    if (y <= 0) return 0.;
    else if (y >= 1) return 0.;
    else return 3 * y * y * (1 - y) * (1 - y) / (den * den);
}

double d2Kijdy2(double y) {
    double den = 3 * y * y - 3 * y + 1;

    if (y <= 0) return 0.;
    else if (y >= 1) return 0.;
    else return (12 * y * y * y - 18 * y * y + 6 * y) / (den * den * den);
}



