#ifndef INFOBLADE_H
#define INFOBLADE_H

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iterator>
#include <ratio>
#include <chrono>
#include <string>
#include <thread>   
#include <vector>
#include <complex>
#include <random>
#include <algorithm>   
#include <cstdio>

#include "../../../include/solver1.h"    
#include "../../../include/blockMeshGenerator.h"
#include "../../../include/sqlite3.h"
#include "../../../include/pathFindeer.h"

#define e 2.718281828459045
//#define PI 3.14159265
#define RadToDegree 57.29577951
#define Cp 1004.5                               
#define gamma 1.4

template<typename T>
using dVec = std::vector<std::vector<T>>;

template<typename T>
using sVec = std::vector<T>;


constexpr double PI = 3.14159265;

namespace infoBlade
{
    //output file stream
    extern std::ofstream fileOut;
    extern std::ofstream fileOut2;
    extern std::ofstream fileOut3;
    extern std::ofstream fileOut4;
    extern std::ofstream compShape;
    extern std::ofstream batchAnalysis[11];

    //design parameters
    extern double initHub;
    extern double VX[3];
    extern dVec<double> tip_radi;
    extern dVec<double> hub_radi;
    extern dVec<double> mean_radi;
    extern double omega1;
    extern double omega2;              
    extern dVec<double> v;
    extern sVec<double> work;
    extern dVec<double> diffusion;
    extern dVec<double> numBlades;
    extern int lowSize;
    extern int highSize;
    extern int totalSize;

    //aerodynamics parameters

    //Mach[i][0] is for velocity 2 relative
    //Mach[i][1] is for velocity 3 absolute 
    extern dVec<std::vector<double>> Mach;
    extern dVec<std::vector<double>> alpha;
    extern dVec<std::vector<double>> beta;
    extern dVec<double> meanAlpha;
    extern dVec<double> meanBeta;
    extern double Mach0;
    extern double resolution;
    extern sVec<double> PR;
    extern sVec<double> R;
    extern dVec<double> Ry;
    extern dVec<double> Area;
    extern dVec<double> chord;
    extern dVec<double> maxCam;       
    extern dVec<double> maxCamPos;
    extern dVec<std::vector<double>> liftCoefficient;
    extern dVec<std::vector<double>> rotateAngle;
    extern dVec<std::vector<double>> incidenceAngle;
    extern dVec<std::vector<double>> AoA;
    //dummy variables
    extern double dummyMaxCam, dummyMaxCamPos, dummyChord, dummyLiftCoefficient;
    //work-done factor
    static double WDF[18] = { 0.982, 0.952, 0.929, 0.910, 0.895, 0.882, 0.875, 0.868, 0.863, 0.860, 0.857, 0.855, 0.853, 0.851, 0.850, 0.849, 0.848, 0.847 };
    extern dVec<std::vector<double>> lossCoefficient;
    extern dVec<std::vector<double>> pressureLoss;
    extern dVec<std::vector<double>> dischargeAngle;
    static double aValues[5] = { 0.2969, -0.1260, -0.3516, 0.2843, -0.1015 } ;
    extern dVec<double> solidity;


    //thermodynamics parameters
    extern double T1;                          
    extern double rho1;
    extern double P1;
    extern dVec<double> Temperature;
    extern dVec<double> TemperatureStag;
    extern dVec<double> Pressure;
    extern dVec<double> PressureStag;
    extern dVec<double> rho;
    extern sVec<double> psi;
    extern sVec<double> phi;
    extern sVec<double> a;
    extern sVec<double> b;
    extern sVec<double> Wr;
    extern sVec<double> Ws;
    extern dVec<double> efficiency;
    
    //misc
    extern double h, g;        

    bool dataBaseSetUp();
    bool initConditionSetUp();    
    void storeInThermoDatabase(sqlite3* db, std::string variablesName, double value, int stage);
    void storeInAeroDatabase(sqlite3* db, std::string variablesName, double value, int stage);
    void storeInDesignDatabase(sqlite3* db, std::string variablesName, double value, int stage);


};

#endif