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


class infoBlade
{
    public : 
    //output file stream
    static std::ofstream fileOut;
    static std::ofstream fileOut2;
    static std::ofstream fileOut3;
    static std::ofstream fileOut4;
    static std::ofstream compShape;
    static std::ofstream batchAnalysis[11];

    //design parameters
    static double VX[3];
    static dVec<double> tip_radi;
    static dVec<double> hub_radi;
    static dVec<double> mean_radi;
    static double omega1;
    static double omega2;              
    static dVec<double> v;
    static sVec<double> work;
    static dVec<double> diffusion;
    static dVec<double> numBlades;
    static int lowSize;
    static int highSize;
    static int totalSize;

    //aerodynamics parameters

    //Mach[i][0] is for velocity 2 relative
    //Mach[i][1] is for velocity 3 absolute 
    static dVec<std::vector<double>> Mach;
    static dVec<std::vector<double>> alpha;
    static dVec<std::vector<double>> beta;
    static dVec<double> meanAlpha;
    static dVec<double> meanBeta;
    static double Mach0;
    static double resolution;
    static sVec<double> PR;
    static sVec<double> R;
    static dVec<double> Area;
    static dVec<double> chord;
    static dVec<double> maxCam;       
    static dVec<double> maxCamPos;
    static dVec<std::vector<double>> liftCoefficient;
    static dVec<std::vector<double>> rotateAngle;
    static dVec<std::vector<double>> incidenceAngle;
    static dVec<std::vector<double>> AoA;
    //dummy variables
    static double dummyMaxCam, dummyMaxCamPos, dummyChord, dummyLiftCoefficient;
    //work-done factor
    double WDF[18] = { 0.982, 0.952, 0.929, 0.910, 0.895, 0.882, 0.875, 0.868, 0.863, 0.860, 0.857, 0.855, 0.853, 0.851, 0.850, 0.849, 0.848, 0.847 };
    static dVec<std::vector<double>> lossCoefficient;
    static dVec<std::vector<double>> pressureLoss;
    static dVec<std::vector<double>> dischargeAngle;
    double aValues[5] = { 0.2969, -0.1260, -0.3516, 0.2843, -0.1015 } ;
    static dVec<double> solidity;


    //thermodynamics parameters
    static double T1;                          
    static double rho1;
    static double P1;
    static dVec<double> Temperature;
    static dVec<double> TemperatureStag;
    static dVec<double> Pressure;
    static dVec<double> PressureStag;
    static dVec<double> rho;
    static sVec<double> psi;
    static sVec<double> phi;
    static sVec<double> a;
    static sVec<double> b;
    static sVec<double> Wr;
    static sVec<double> Ws;
    static dVec<double> efficiency;
    
    //misc
    static double h, g;        

    bool dataBaseSetUp();
    bool initConditionSetUp();                  

};

#endif