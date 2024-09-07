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

    //output file stream
    std::ofstream fileOut;
    std::ofstream fileOut2;
    std::ofstream fileOut3;
    std::ofstream fileOut4;
    std::ofstream compShape;
    std::ofstream batchAnalysis[11];

    //design parameters
    double VX[3];
    dVec<double> tip_radi;
    dVec<double> hub_radi;
    dVec<double> mean_radi;
    double omega1;
    double omega2;              
    dVec<double> v;
    sVec<double> work;
    dVec<double> diffusion;
    dVec<double> numBlades;
    int lowSize;
    int highSize;
    int totalSize;

    //aerodynamics parameters

    //Mach[i][0] is for velocity 2 relative
    //Mach[i][1] is for velocity 3 absolute 
    dVec<std::vector<double>> Mach;
    dVec<std::vector<double>> alpha;
    dVec<std::vector<double>> beta;
    dVec<double> meanAlpha;
    dVec<double> meanBeta;
    double Mach0;
    double resolution;
    sVec<double> PR;
    sVec<double> R;
    dVec<double> Area;
    dVec<double> chord;
    dVec<double> maxCam;       
    dVec<double> maxCamPos;
    dVec<std::vector<double>> liftCoefficient;
    dVec<std::vector<double>> rotateAngle;
    dVec<std::vector<double>> incidenceAngle;
    dVec<std::vector<double>> AoA;
    //dummy variables
    double dummyMaxCam, dummyMaxCamPos, dummyChord, dummyLiftCoefficient;
    //work-done factor
    double WDF[18] = { 0.982, 0.952, 0.929, 0.910, 0.895, 0.882, 0.875, 0.868, 0.863, 0.860, 0.857, 0.855, 0.853, 0.851, 0.850, 0.849, 0.848, 0.847 };
    dVec<std::vector<double>> lossCoefficient;
    dVec<std::vector<double>> pressureLoss;
    dVec<std::vector<double>> dischargeAngle;
    double aValues[5] = { 0.2969, -0.1260, -0.3516, 0.2843, -0.1015 } ;
    dVec<double> solidity;


    //thermodynamics parameters
    double T1;                          
    double rho1;
    double P1;
    dVec<double> Temperature;
    dVec<double> TemperatureStag;
    dVec<double> Pressure;
    dVec<double> PressureStag;
    dVec<double> rho;
    sVec<double> psi;
    sVec<double> phi;
    sVec<double> a;
    sVec<double> b;
    sVec<double> Wr;
    sVec<double> Ws;
    dVec<double> efficiency;
    

    //misc
    double h, g;

    void initConditionSetUp();                  

};