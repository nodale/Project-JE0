#include <cstdio>
#include <filesystem>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iterator>
#include <ratio>
#include <chrono>
#include <thread>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "/home/nodale/Documents/JE0/JE0/include/solver1.h"    
#include "/home/nodale/Documents/JE0/JE0/include/blockMeshGenerator.h"

#define e 2.718281828459045
//#define PI 3.14159265
#define RadToDegree 57.29577951
#define Cp 1004.5                               
#define gamma 1.4

//TODO
//improve the optomiseFlow function, make it pick the best combination

template<typename T>
using dVec = std::vector<std::vector<T>>;

template<typename T>
using sVec = std::vector<T>;


constexpr double PI = 3.14159265;

class Blade : public solver1, public blockMeshGen
{

std::ofstream fileOut;
std::ofstream fileOut2;
std::ofstream fileOut3;
std::ofstream fileOut4;

std::ofstream compShape;

std::ofstream batchAnalysis[11];

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

//thermodynamics initial conditions
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
dVec<std::vector<double>> lossCoefficient;
dVec<std::vector<double>> pressureLoss;
dVec<std::vector<double>> dischargeAngle;
double aValues[5] = { 0.2969, -0.1260, -0.3516, 0.2843, -0.1015 } ;

dVec<double> solidity;


double h, g;

public:
Blade(dVec<double> tipP, double hub, double w1, double w2, double res, double vx1, double vx2, sVec<double> deltaP, double Temp, double P, sVec<double> initial_alpha1, sVec<double> Reaction, dVec<double> c, int lpSize, int hpSize)
: omega1 { w1 }
, omega2 { w2 }
, resolution { res }
, VX { vx1, vx2, vx1 }
, T1 { Temp }
, P1 { P }
, lowSize( lpSize )     
, highSize( hpSize )
{
    fileOut.open("output/misc/out.dat");
    fileOut2.open("output/misc/out2.dat");
    fileOut3.open("output/misc/out3.dat");  
    fileOut4.open("output/misc/out4.dat");
    compShape.open("output/misc/shape.dat");
    char camFilename[] = "output/misc/camberline.dat";
    solver1::initFile(camFilename);

    totalSize = lpSize + hpSize;

    tip_radi.resize(totalSize);
    hub_radi.resize(totalSize);
    mean_radi.resize(totalSize);
              
    v.resize(totalSize);
    work.resize(totalSize);
    diffusion.resize(totalSize);
    numBlades.resize(totalSize);

    Mach.resize(totalSize);

    alpha.resize(totalSize + 1);
    beta.resize(totalSize);
    meanAlpha.resize(lpSize + hpSize + 1);
    meanBeta.resize(totalSize);
    PR.resize(totalSize);
    R.resize(totalSize);
    Area.resize(totalSize + 1);
    chord.resize(totalSize);
    maxCam.resize(totalSize);       
    maxCamPos.resize(totalSize);
    liftCoefficient.resize(totalSize);
    rotateAngle.resize(totalSize);
    incidenceAngle.resize(totalSize);
    AoA.resize(totalSize);

    Temperature.resize(totalSize);
    TemperatureStag.resize(totalSize);
    Pressure.resize(totalSize);
    PressureStag.resize(totalSize);
    rho.resize(totalSize);
    psi.resize(totalSize);
    phi.resize(totalSize);
    a.resize(totalSize);
    b.resize(totalSize);
    Wr.resize(totalSize);
    Ws.resize(totalSize);
    efficiency.resize(totalSize);          
    lossCoefficient.resize(totalSize);
    pressureLoss.resize(totalSize);
    dischargeAngle.resize(totalSize);

    solidity.resize(totalSize);

    for(int i = 0; i < totalSize; i++)
    {
        tip_radi[i].resize(3);
        tip_radi[i].resize(3);
        hub_radi[i].resize(3);
        mean_radi[i].resize(3);
        v[i].resize(3);
        diffusion[i].resize(2);
        numBlades[i].resize(2);

        Mach[i].resize(2);

        alpha[i].resize(2);
        beta[i].resize(2);
        meanAlpha[i].resize(2);
        meanBeta[i].resize(2);
        Area[i].resize(3);
        chord[i].resize(2);
        maxCam[i].resize(2);       
        maxCamPos[i].resize(2);
        liftCoefficient[i].resize(2);
        rotateAngle[i].resize(2);
        incidenceAngle[i].resize(2);
        AoA[i].resize(2);
        
        Temperature[i].resize(3);
        TemperatureStag[i].resize(3);
        Pressure[i].resize(3);
        PressureStag[i].resize(3);
        rho[i].resize(3);
        efficiency[i].resize(2);       
        lossCoefficient[i].resize(2);
        pressureLoss[i].resize(2);
        dischargeAngle[i].resize(2);
        
        solidity[i].resize(2);
    }

    alpha[alpha.size() - 1].resize(2);

    std::string tempName;   
    for(int i = 0; i < totalSize; i++)
    {
        tip_radi[i] = tipP[i];

        PR[i] = deltaP[i];  
        meanAlpha[i][0] = initial_alpha1[i];    
        R[i] = Reaction[i];

        for(int j = 0; j < 2; j++)
        {
            diffusion[i][j] = 0.5;  

            chord[i][0] = c[i][0];      
            chord[i][1] = c[i][1];

            alpha[i][j].resize(resolution);
            beta[i][j].resize(resolution);                                         
            lossCoefficient[i][j].resize(resolution);
            liftCoefficient[i][j].resize(resolution);      
            Mach[i][j].resize(resolution);
            pressureLoss[i][j].resize(resolution);         
            dischargeAngle[i][j].resize(resolution);
            rotateAngle[i][j].resize(resolution);
            incidenceAngle[i][j].resize(resolution);
            AoA[i][j].resize(resolution);
            //this doesn't work for some reasons
            //efficiency[i][j] = 0.9;s
            
        }

        //opening a file for each batchAnalysis ofstream
        tempName = "batchData/" + std::to_string(i) + ".dat";
        batchAnalysis[i].open( tempName);
    }

    alpha[alpha.size() - 1][0].resize(resolution);
    alpha[alpha.size() - 1][1].resize(resolution);
    
    for(int i = 0; i <= resolution; i++)
    {
        alpha[alpha.size() - 1][0][i] = initial_alpha1[alpha.size() - 1];
        alpha[alpha.size() - 1][1][i] = initial_alpha1[alpha.size() - 1];

    }

    Area[Area.size()-1].resize(2);
    meanAlpha[meanAlpha.size()-1].resize(2);
    meanAlpha[meanAlpha.size()-1][0] = initial_alpha1[meanAlpha.size()-1];

    hub_radi[0][0] = hub;
    mean_radi[0][0] = 0.5 * ( tip_radi[0][0] + hub_radi[0][0] );            

    Temperature[0][0] = T1;
    Pressure[0][0] = P1;

    //initialising data for stage 1
    //Stagnation temperature 1 need mean alpha
    TemperatureStag[0][0] = Temperature[0][0] + 0.5 * pow( VX[0] / cos( meanAlpha[0][0] / RadToDegree ) , 2 ) / Cp ;
    //phi
    phi[0] =  VX[0] / ( omega1 * mean_radi[0][0] );
    //beta 1
    meanBeta[0][0] = atan( tan( meanAlpha[0][0] / RadToDegree ) - 1 / phi[0] ) * RadToDegree;
    //Temperature 3
    Temperature[0][2] = Temperature[0][0] * pow( ( PR[0] ) , ( gamma - 1 ) / gamma );
    //Stagnation temperature 3 need mean alpha
    TemperatureStag[0][2] = Temperature[0][2] + 0.5 * pow( VX[0] / cos( meanAlpha[1][0] / RadToDegree ) , 2 ) / Cp;
    //Stagnation pressure  1
    PressureStag[0][0] = Pressure[0][0] * pow( ( TemperatureStag[0][0] / Temperature[0][0] ) , gamma / ( gamma - 1 ) );
    //Stagnation pressure 3
    PressureStag[0][2] = PressureStag[0][0] * PR[0];
    //Temperature stagnation 2
    work[0] = Cp * ( TemperatureStag[0][2] - TemperatureStag[0][0] );
    TemperatureStag[0][1] = TemperatureStag[0][0] + ( R[0] * work[0] / Cp );
    //Pressure 2
    Pressure[0][1] = Pressure[0][0] * pow( ( TemperatureStag[0][1] / TemperatureStag[0][0] ) , gamma / ( gamma - 1 ) );
    //Stagnation pressure 2
    PressureStag[0][1] = PressureStag[0][0] * pow( ( TemperatureStag[0][1] / TemperatureStag[0][0] ) , gamma / ( gamma - 1 ) );
    //psi
    psi[0] = work[0] / pow( omega1 * mean_radi[0][0] , 2 );
    //alpha 2
    meanAlpha[0][1] = atan2( psi[0] , phi[0] * WDF[0] ) * RadToDegree;
     //beta 2
    meanBeta[0][1] = atan( tan( meanAlpha[0][1] / RadToDegree ) - 1 / phi[0] ) * RadToDegree;
    //Temperature 2
    Temperature[0][1] = TemperatureStag[0][1] - 0.5 * pow(  VX[0] / cos( meanAlpha[0][1] / RadToDegree ) , 2 ) / Cp;
    //Pressure 2
    Pressure[0][1] = Pressure[0][0] * pow( ( Temperature[0][1] / Temperature[0][0] ) , gamma / ( gamma - 1 ) );
    //Pressure 3 need mean alpha
    Pressure[0][2] = PressureStag[0][2] * pow( ( TemperatureStag[0][2] / Temperature[0][2] ) , -gamma / ( gamma - 1 ) );
   //Area stag 1, station 1
    Area[0][0] = PI * (  pow( tip_radi[0][0] , 2 ) - pow( hub_radi[0][0] , 2 ) );

    //rho 1,2,3
    for(int i = 0; i < 3; i++)
    {
    rho[0][i] = 1000 * Pressure[0][i] / ( 287 * Temperature[0][i] );
    }

    for(int i = 1; i < 3; i++)
    {
    Area[0][i] = Area[0][i-1] * ( rho[0][i-1] / rho[0][i] );
    }

    Area[1][0] = Area[0][2];

    for(int i = 0; i < 3; i++)
    {
    hub_radi[0][i] = sqrt( -( Area[0][i] / PI ) + pow( tip_radi[0][i] , 2 )  ); 
    }

    for(int i = 1; i < 3; i++)
    {
    mean_radi[0][i] = 0.5 * ( tip_radi[0][0] + hub_radi[0][0] );
    }

    //solidity[0] = ( 1.5 * psi[0] ) / ( 1.55 * phi[0] - psi[0] );
    // solidity[0][0] = 0.5 * ( ( tan(meanBeta[0][0]/RadToDegree) - tan(meanBeta[0][1]/RadToDegree) ) ) / ( ( diffusion[0][0] - ( 1.0 - ( cos(meanBeta[0][0]/RadToDegree) / cos( meanBeta[0][1]/RadToDegree ) ) ) ) / cos(meanBeta[0][0]/RadToDegree) ); 
    // solidity[0][1] = 0.5 * ( ( tan(meanAlpha[0][1]/RadToDegree) - tan(meanAlpha[1][0]/RadToDegree) ) ) / ( ( diffusion[0][1] - ( 1.0 - ( cos(meanAlpha[1][0]/RadToDegree) / cos( meanAlpha[0][1]/RadToDegree ) ) ) ) / cos(meanAlpha[0][1]/RadToDegree) ); 
    solidity[0][0] = fabs( ( 1.5 * psi[0] ) / ( 1.55 * phi[0] - psi[0] ) );
    solidity[0][1] = fabs( ( 1.5 * psi[0] ) / ( 1.55 * phi[0] - psi[0] ) );


    numBlades[0][0] = 2 * PI * mean_radi[0][0] / ( chord[0][0] / solidity[0][0] );  
    numBlades[0][1] = 2 * PI * mean_radi[0][1] / ( chord[0][1] / solidity[0][1] );

    chord[0][0] = solidity[0][0] * 2 * PI * mean_radi[0][0] / numBlades[0][0];
    chord[0][1] = solidity[0][1] * 2 * PI * mean_radi[0][1] / numBlades[0][1];

    double vu1_r = VX[0] * tan( meanAlpha[0][0] / RadToDegree );

    //1st approach
    //double vu2_r = VX[0] * tan( meanAlpha[0][1] / RadToDegree );

    //2nd approach
    double vu2_r = ( work[0] + omega1 * vu1_r * mean_radi[0][0] ) / ( omega1 * mean_radi[0][0] );
    a[0] = 0.5 * ( vu1_r + vu2_r );
    b[0] = 0.5 * ( vu2_r - vu1_r );

    //std::cout << Pressure[0][0] << " " << Pressure[0][1] << " " << Pressure[0][2] << std::endl;

}

//to initialise and calculate all thermodynamics and other variables from stage 2 to 11
void init()
{
    for(int i = 1; i < totalSize; i++)
    {
    //copying data from station 3 of the previous stage to station 1 of current stage
    TemperatureStag[i][0] = TemperatureStag[i-1][2];
    PressureStag[i][0] = PressureStag[i-1][2];
    Temperature[i][0] = Temperature[i-1][2];
    Pressure[i][0] = Pressure[i-1][2];
    rho[i][0] = rho[i-1][0];    
    mean_radi[i][0] = mean_radi[i-1][2];
    hub_radi[i][0] = hub_radi[i-1][2];
    
    //phi
    if( i < lowSize )
    {
    phi[i] =  VX[0] / ( omega1 * mean_radi[i][0] );
    }
    if( i >= lowSize )
    {
    phi[i] =  VX[0] / ( omega2 * mean_radi[i][0] );
    }

    //beta 1
    meanBeta[i][0] = atan( tan( meanAlpha[i][0] / RadToDegree ) - 1 / phi[i] ) * RadToDegree;
    //Temperature 3
    Temperature[i][2] = Temperature[i][0] * pow( ( PR[i] ) , ( gamma - 1 ) / gamma );
    //Stagnation temperature 3, needs mean alpha
    TemperatureStag[i][2] = Temperature[i][2] + 0.5 * pow( VX[0] / cos( meanAlpha[i+1][0] / RadToDegree ) , 2 ) / Cp;
    //Stagnation pressure 3
    PressureStag[i][2] = PressureStag[i][0] * PR[i];
    //Temperature stagnation 2
    work[i] = Cp * ( TemperatureStag[i][2] - TemperatureStag[i][0] );
    TemperatureStag[i][1] = TemperatureStag[i][0] + ( R[i] * work[i] / Cp );
    //Pressure 2
    Pressure[i][1] = Pressure[i][0] * pow( ( TemperatureStag[i][1] / TemperatureStag[i][0] ) , gamma / ( gamma - 1 ) );
    //Stagnation pressure 2
    PressureStag[i][1] = PressureStag[i][0] * pow( ( TemperatureStag[i][1] / TemperatureStag[i][0] ) , gamma / ( gamma - 1 ) );
    //psi
    if( i < lowSize ) 
    {
    psi[i] = work[i] / pow( omega1 * mean_radi[i][0] , 2 );
    }
    if( i >= lowSize )
    {
    psi[i] = work[i] / pow( omega2 * mean_radi[i][0] , 2 );
    }
    //alpha 2
    meanAlpha[i][1] = atan2( psi[i] , phi[i] ) * RadToDegree;
     //beta 2
     
    meanBeta[i][1] = atan( tan( meanAlpha[i][1] / RadToDegree ) - 1 / phi[i] ) * RadToDegree;
    //Temperature 2
    Temperature[i][1] = TemperatureStag[i][1] - 0.5 * pow(  VX[0] / cos( meanAlpha[i][1] / RadToDegree ) , 2 ) / Cp;
    //Pressure 2
    Pressure[i][1] = Pressure[i][0] * pow( ( Temperature[i][1] / Temperature[i][0] ) , gamma / ( gamma - 1 ) );
    //Pressure 3, needs mean alpha
    Pressure[i][2] = PressureStag[i][2] * pow( ( TemperatureStag[i][2] / Temperature[i][2] ) , -gamma / ( gamma - 1 ) );
   //Area station 1
    Area[i][0] = Area[i-1][2];

    //rho 1,2,3
    for(int j = 0; j < 3; j++)
    {
    rho[i][j] = 1000 * Pressure[i][j] / ( 287 * Temperature[i][j] );
    }

    for(int j = 1; j < 3; j++)
    {
    Area[i][j] = Area[i][j-1] * ( rho[i][j-1] / rho[i][j] );
    }

    Area[i+1][0] = Area[i][2];

    for(int j = 0; j < 3; j++)
    {
    hub_radi[i][j] = pow( -( Area[i][j] / PI ) + pow( tip_radi[i][j] , 2 ) , 0.5  ); 
    }

    for(int j = 1; j < 3; j++)
    {
    mean_radi[i][j] = 0.5 * ( tip_radi[i][j] + hub_radi[i][j] );                                
    }

    //solidity[i] = fabs( ( 1.5 * psi[i] ) / ( 1.55 * phi[i] - psi[i] ) );
    // solidity[i][0] = 0.5 * ( ( tan(meanBeta[i][0]/RadToDegree) - tan(meanBeta[i][1]/RadToDegree) ) ) / ( ( diffusion[i][0] - ( 1.0 - ( cos(meanBeta[i][0]/RadToDegree) / cos( meanBeta[i][1]/RadToDegree ) ) ) ) / cos(meanBeta[i][0]/RadToDegree) );   
    // solidity[i][1] = 0.5 * ( ( tan(meanAlpha[i][1]/RadToDegree) - tan(meanAlpha[i+1][0]/RadToDegree) ) ) / ( ( diffusion[i][1] - ( 1.0 - ( cos(meanAlpha[i+1][0]/RadToDegree) / cos( meanAlpha[i][1]/RadToDegree ) ) ) ) / cos(meanAlpha[i][1]/RadToDegree) ); 

    solidity[i][0] = fabs( ( 1.5 * psi[i] ) / ( 1.55 * phi[i] - psi[i] ) );
    solidity[i][1] = fabs( ( 1.5 * psi[i] ) / ( 1.55 * phi[i] - psi[i] ) );


    double vu1_r = VX[0] * tan( meanAlpha[i][0] / RadToDegree );
    // approach 1
    //double vu2_r = VX[0] * tan( meanAlpha[i][1] / RadToDegree );
    
    //approach 2
    double vu2_r;
    if(i < lowSize)
    {
        vu2_r = ( work[i] + omega1 * vu1_r * mean_radi[i][0] ) / ( omega1 * mean_radi[i][0] );
    }
    if(i >= lowSize)
    {
        vu2_r = ( work[i] + omega2 * vu1_r * mean_radi[i][0] ) / ( omega2 * mean_radi[i][0] );
    }

    numBlades[i][0] = 2 * PI * mean_radi[i][0] / ( chord[i][0] / solidity[i][0] );
    numBlades[i][1] = 2 * PI * mean_radi[i][1] / ( chord[i][1] / solidity[i][1] );

    chord[i][0] = solidity[i][0] * 2 * PI * mean_radi[i][0] / numBlades[i][0];
    chord[i][1] = solidity[i][1] * 2 * PI * mean_radi[i][1] / numBlades[i][1];

    //std::cout << TemperatureStag[i][0] << " " << PressureStag[i][0]  << " " << PressureStag[i][2] << std::endl;

    a[i] = 0.5 * ( vu1_r + vu2_r );
    b[i] = 0.5 * ( vu2_r - vu1_r );

    }


    //for some reasons chord[0][0] changes here
    chord[0][0] = chord[1][0];
    //std::cout << chord[0][0] << " " << chord[0][1] << std::endl;
    //std::cout << Pressure[10][2] << std::endl;        
}

void getFlowPaths(int i)
{
    
    double dr1, dr2, radius1, radius2, y;
    double tempPhi, tempVu1, tempVu2;
    dr1 = ( tip_radi[i][0] - hub_radi[i][0] ) / resolution ;
    dr2 = ( tip_radi[i][1] - hub_radi[i][1] ) / resolution ;

    for(int r = 0; r <= resolution; r++)
    {
        radius1 = hub_radi[i][0] + r * dr1;
        radius2 = hub_radi[i][1] + r * dr2;
        
        y = radius1 / mean_radi[i][0]; 

        tempVu1 = (a[i] - b[i] ) /y;
        tempVu2 = (a[i] + b[i] ) /y;

        //approach from energy method, with angle 0 at the leading edge
        // tempVu1 = 0.0;
        // if(i < 3)
        // {
        // tempVu2 = work[i] / (omega1 * radius2);
        // }
        // if(i >= 3)
        // {
        // tempVu2 = work[i] / (omega2 * radius2);
        // }
        
        //TODO get phi for every radius
        if( i < lowSize)
        {
        tempPhi = VX[0] / ( omega1 * radius1 );
        }
        if(i >= lowSize)
        {
        tempPhi = VX[0] / ( omega2 * radius1 );  
        }

        alpha[i][0][r] = atan( tempVu1 / VX[0] ) * RadToDegree;
        alpha[i][1][r] = atan( tempVu2 / VX[0] ) * RadToDegree;

        beta[i][0][r] = atan( tan( alpha[i][0][r] / RadToDegree ) - 1 / tempPhi ) * RadToDegree;
        beta[i][1][r] = atan( tan( alpha[i][1][r] / RadToDegree ) - 1 / tempPhi ) * RadToDegree;

        //alpha[11][0][r] = meanAlpha[11][0];
        //printOut(fileOut, r, alpha[i][0][r]); 
        //printOut(fileOut2, r, alpha[i][1][r]);
        //printOut(fileOut3, r, beta[i][0][r]);
        //printOut(fileOut4, r, beta[i][1][r]);       
    }
}

//j is either one(rotor) or two(stator)
void getCamberline(int i, int j, int r)
{
    double dummyR = 0.0;
    //distance to the maximum camber, dimesionless
    maxCamPos[i][j] = 0.5 * chord[i][j];
    double theta;

    if(j == 0)
    {
    theta = beta[i][0][r] - beta[i][1][r]; 
    }
    if(j == 1)
    {
    theta = alpha[i][0][r] - alpha[i][1][r]; 
    }

    maxCam[i][j] = chord[i][j] / ( 4 * tan( theta / RadToDegree ) ) * ( sqrt( fabs( 1 + pow( 4 * ( tan( theta / RadToDegree ) ) , 2 ) * ( maxCamPos[i][j] / chord[i][j] - pow( maxCamPos[i][j] / chord[i][j] - 3.0 / 16.0 , 2 ) ) ) ) - 1 );

    //initialising dummy variables
    dummyMaxCam = maxCam[i][j];
    dummyMaxCamPos = maxCamPos[i][j];
    dummyChord = chord[i][j]; 


    //std::cout << theta << " " << dummyMaxCam << " " << dummyMaxCamPos << " " << dummyChord << std::endl;

    for(double dt = 0; dt < dummyChord;)
    {
        dummyR = func( dt, dummyR );
        dt += dummyChord/resolution;
        file_out << dt << " " << dummyR << "\n";
    }

    //solver1::rungeKuttam(dummyChord, 0.0, 0.01 / resolution, 0.0, dummyChord, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);

}

void analysisCamber(int j)
{
    for(int i = 0; i < totalSize; i++)
    {
    //maxCamPos[i][j] = 0.5 * chord[i][j];
    double theta, yValue;

    //double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;

    for(int r = 0; r <= resolution; r++)
    {
    // // //yValue = ( hub_radi[i][j] + dr * r ) / mean_radi[i][j];
    if(j == 0)
    {
    theta = cos (beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    }
    if(j == 1)
    {
    theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    }

    // if(j == 0)s
    // {
    // theta = beta[i][0][r] - beta[i][1][r]; 
    // }
    // if(j == 1)
    // {
    // theta = alpha[i][1][r] - alpha[i+1][0][r]; 
    // }

    //maxCam[i][j] = chord[i][j] / ( 4 * tan( theta / RadToDegree ) ) * ( sqrt( fabs( 1 + pow( 4 * ( tan( theta / RadToDegree ) ) , 2 ) * ( maxCamPos[i][j] / chord[i][j] - pow( maxCamPos[i][j] / chord[i][j] - 3.0 / 16.0 , 2 ) ) ) ) - 1 );

    // if(j == 1)
    // {
    // batchAnalysis[i] << r << " " << alpha[i][j][r] << "\n";
    // }
    // if(j == 0)
    // {
    // batchAnalysis[i] << r << " " << beta[i][j][r] << "\n";
    // }

    //get loss coefficients, not extremely accurate but it's alright
    lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i][1] / cos( alpha[i+1][0][r] / RadToDegree );

    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    //pressure loss is in Kpascal
    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

    batchAnalysis[i] << theta << " " <<  liftCoefficient[i][j][r] << "\n";

    //efficiency analysis
    //double dummyEfficiency = 1 - ( lossCoefficient[i][0][r] / pow ( cos( alpha[i+1][0][r] / RadToDegree ), 2 ) + lossCoefficient[i][1][r] / pow ( cos( beta[i][0][r] / RadToDegree ), 2 ) ) * pow( phi[i] , 2 ) / ( 2 * psi[i] );
    //batchAnalysis[i] << r << " " <<  dummyEfficiency << "\n";
    }
    }
    //solver1::rungeKuttam(dummyChord, 0.0, 0.01 / resolution, 0.0, dummyChord, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);

    //drawShape();
}

void getLeadingTrailingAngles(int j)
{
    double dt, discharge1, discharge2, dummyDischargeAngle, dummyR;
    for(int i = 0; i < totalSize; i++)
    {
        dummyMaxCamPos = 0.5 * chord[i][j];
        dummyChord = chord[i][j];

        for(int r = 0; r <= resolution; r++)
        {
            dt = dummyChord / resolution;
            dummyLiftCoefficient = liftCoefficient[i][j][r];
            discharge1 = func( dt * ( resolution - 2 ) , discharge1 );
            discharge2 = func( dt * ( resolution - 1 ) , discharge2 );

            dischargeAngle[i][j][r] = atan2(  discharge2 - discharge1  , dt ) * RadToDegree;
            
            if( j == 0)
            {
            rotateAngle[i][j][r] = beta[i][1][r] - dischargeAngle[i][j][r];
            incidenceAngle[i][j][r] = beta[i][0][r] - rotateAngle[i][j][r];
            //batchAnalysis[i] << beta[i][0][r] << " " <<  rotateAngle[i][j][r] << "\n";
            }
            if( j == 1)
            {
            rotateAngle[i][j][r] = alpha[i+1][0][r] - dischargeAngle[i][j][r];
            incidenceAngle[i][j][r] = alpha[i][1][r] - rotateAngle[i][j][r];
            //batchAnalysis[i] << r << " " <<  incidenceAngle[i][j][r] << "\n";
            }
            //batchAnalysis[i] << r << " " <<  AoA[i][j][r] << "\n";
        }
    }
}

void getBladeAngles(int i, int j)
{
    double dr1, dr2, radius1, radius2, E, alphaPoint;
    double theta, k1, k2, dummyK1, dummyK2, delta, dIncidence, incidencePoint, deltaPoint, Kti, Ktdelta;
    double q, delta0Point, m, m10, tempX, tempB;

    dr1 = ( tip_radi[i][0] - hub_radi[i][0] ) / resolution ;
    dr2 = ( tip_radi[i][1] - hub_radi[i][1] ) / resolution ;

    //chosen incidence angle
    double incidence;
    double stagger;
    double Ksh = 1.0;
    double tb = 0.1 * chord[i][0];
    //the distance to maximum thickness, might need to change it later
    double distToMax = 0.5;
    bool stop = 0;
    int correctionVar = 1;
    //everything is positive, need to readjust them again later
    for(int r = 0; r <= resolution; r++)
    {
        //stagger angle, currently it is the difference between inlet and outlet angle
        incidence = incidenceAngle[i][j][r];
        if(j == 0)
        {
            stagger = fabs( beta[i][1][r] - beta[i][0][r] );
            k1 = fabs(beta[i][0][r]);
            k2 = fabs(beta[i][1][r]);
        }
        if(j == 1)
        {
            stagger = fabs( alpha[i+1][0][r] - alpha[i][1][r] );
            k1 = fabs(alpha[i][0][r]);
            k2 = fabs(alpha[i][1][r]);
        }
        stop = 0;
        radius1 = hub_radi[i][0] + r * dr1;
        radius2 = hub_radi[i][1] + r * dr2;
        

        //need a loop for the following, use the stop boolean and if statement commented below
        for(int count = 0; count <= resolution; count++)
        {
            dummyK1 = k1;
            dummyK2 = k2;
            theta = k1 - k2; 

            if(theta < 0)
            {
            correctionVar = -1;
            }
            if(theta >= 0)
            {
            correctionVar = 1;
            }

            E = 0.65 - 0.002 * theta;
            q = 0.28 / ( 0.1 + pow( tb / chord[i][0] , 0.3 ) );
            Kti = pow( 10 * tb / chord[i][0] , q );
            alphaPoint = correctionVar * ( 3.6 * Ksh * Kti + 0.3532 * theta * pow( distToMax / chord[i][0] , 0.25 ) ) * pow( correctionVar * theta , E );
            incidencePoint = alphaPoint + stagger - k1;

            if(j == 0)
            {
            delta0Point = 0.01 * solidity[i][0] * beta[i][0][r] + ( 0.74 * pow( solidity[i][0] , 1.9 ) + 3 * solidity[i][0] ) * pow( fabs( beta[i][0][r] ) / 90 , 1.67 + 1.09 * solidity[i][0] );
            }
            if(j == 1)
            {
            delta0Point = 0.01 * solidity[i][1] * beta[i][0][r] + ( 0.74 * pow( solidity[i][1] , 1.9 ) + 3 * solidity[i][1] ) * pow( fabs( alpha[i][1][r] ) / 90 , 1.67 + 1.09 * solidity[i][1] );
            }

            tempX = beta[i][0][r] / 100;
            m10 = 0.17 - 0.0333 * tempX + 0.333 * pow( tempX , 2 );
            tempB = 0.9625 - 0.17 * tempX - 0.85 * pow( tempX , 3 );
            m = m10 / pow( solidity[i][j] , tempB );
            Ktdelta = 6.25 * ( tb / chord[i][0] ) + 37.5 * pow( tb / chord[i][0] , 2 );
            deltaPoint = Ksh * Ktdelta * delta0Point + m * theta;

            dIncidence = incidence - incidencePoint;
            //TODO fix incidencePoint
            if(j == 0)
            {
            k1 = fabs(beta[i][0][r]) - incidencePoint - dIncidence;
            k2 = fabs(beta[i][1][r]) - deltaPoint;
            }
            if(j == 1)
            {
            k1 = fabs(alpha[i][1][r]) - incidencePoint - dIncidence;
            k2 = fabs(alpha[i+1][0][r]) - deltaPoint;
            }
            
            //std::cout << k2 << " " << k2 << std::endl;

            //convergence test
            //if ( dummyK1 == k1 && k2 == dummyK2 )
            if(count == resolution)
            {
                if(j == 0)
                {
                AoA[i][j][r] = k1 - beta[i][0][r];
                
                //std::cout << "inlet :  "<< k1 << " " << fabs(beta[i][0][r]) << "   outlet :  " << k2 << " " << fabs(beta[i][1][r]) << "\n";
                }
                if(j == 1)
                {
                AoA[i][j][r] = k2 - alpha[i][1][r];

                //std::cout << "inlet :  "<< k1 << " " << fabs(alpha[i][1][r]) << "   outlet :  " << k2 << " " << fabs(alpha[i+1][0][r]) << "\n";
                }
                //batchAnalysis[i] << r << " " <<  AoA[i][j][r] << "\n";
                stop = 1;
            }
        }

    }

    //drawFlowPath();
}

void generateBlade(int i, int j)
{
    double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;
    double radius;
    char filename[] = "blade.stl"; //for blockMesh : blockMeshDict
    blockMeshGen::init(filename);

    for(int r = 0; r <= resolution; r++)
    {
        dummyLiftCoefficient = liftCoefficient[i][j][r];
    
        radius = dr * r + hub_radi[i][j];

        //distance to the maximum camber, dimesionless
        dummyChord = chord[i][j]; 
        maxCamPos[i][j] = 0.50 * dummyChord;
        double theta;

        if(j == 0)
        {
        theta = beta[i][0][r] - beta[i][1][r]; 
        }
        if(j == 1)
        {
        theta = alpha[i][1][r] - alpha[i+1][0][r]; 
        }

        maxCam[i][j] = dummyChord / ( 4 * tan( theta / RadToDegree ) ) * ( sqrt(  1 + pow( 4 * ( tan( theta / RadToDegree ) ) , 2 ) * ( maxCamPos[i][j] / dummyChord - pow( maxCamPos[i][j] / dummyChord - 3.0 / 16.0 , 2 ) ) ) - 1 );
        //std::cout << maxCam[i][j] << "\n";
        //initialising dummy variables
        //for(int t = -1; t < 2;)
        //{
        
        dummyMaxCam = maxCam[i][j];
        dummyMaxCamPos = maxCamPos[i][j];
        
        //std::cout << theta << " " << dummyMaxCam << " " << dummyMaxCamPos << " " << dummyChord << std::endl;
        double dummyR = 0.0;
        double dummyDist = 0.0; 
        double dummyT, dummyXU, dummyYU, dummyXL, dummyYL, dummyAngle, dummyPrevX, dummyPrevY;
        int count = 0;
        for(int ra = 0; ra < 2; ra++)
        {
            //UPPER
            for(double dt = 0; dt <= dummyChord;)
            {
                dummyPrevX = dummyDist;
                dummyPrevY = dummyR;
                //calculating camber abscissa for rotated coordinate
                dummyDist = 0.5 * dummyChord - dt;
                //calculating camber ordinate
                dummyR = func( dt, dummyR );
                //calculating thickness yt
                dummyT = 5 * 0.1 * dummyChord * ( sqrt( dummyDist ) * aValues[0]  +  aValues[1] * dummyDist  + aValues[2] * pow( dummyDist , 2 ) + aValues[3] * pow( dummyDist , 3 ) + aValues[4] * pow( dummyDist, 4 ) );
                //getting upper x
                dummyXU = dummyDist - dummyT * sin( atan2( dummyR - dummyPrevY , dummyDist - dummyPrevX ) );
                dummyYU = dummyR + dummyT * cos( atan2( dummyR - dummyPrevY , dummyDist - dummyPrevX ) );

                transformAngle(rotateAngle[i][j][r], dummyXU, dummyYU);//0.5 * dummyChord - dummyDist
                
                blockMeshGen::collectVertices(0.5 * dummyChord - dummyXU, dummyYU, radius);
                
                std::cout << "progress : vertex number " << count << " and radius " << r << " out of " << resolution << std::endl; 
                
                dt += 2 * dummyChord/resolution;
                count += 1;

            }
            // dummyR = 0.0;
            // dummyDist = 0.0;
            // //LOWER
            for(double dt = dummyChord; dt >= 0;)
            {
                //calculating camber abscissa for rotated coordinate
                dummyDist = 0.5 * dummyChord - dt;
                //calculating camber ordinate
                dummyR = func( dt, dummyR );
                //calculating thickness yt
                dummyT = 5 * 0.1 * dummyChord * ( sqrt( dummyDist ) * aValues[0]  +  aValues[1] * dummyDist  + aValues[2] * pow( dummyDist , 2 ) + aValues[3] * pow( dummyDist , 3 ) + aValues[4] * pow( dummyDist, 4 ) );
                
                dummyXL = dummyDist + dummyT * sin( atan2( dummyR - dummyPrevY , dummyDist - dummyPrevX ) );
                dummyYL = dummyR - dummyT * cos( atan2( dummyR - dummyPrevY , dummyDist - dummyPrevX ) );

                transformAngle(rotateAngle[i][j][r], dummyXL, dummyYL);//0.5 * dummyChord - dummyDist
                
                blockMeshGen::collectVertices(0.5 * dummyChord - dummyXL, dummyYL, radius);
                
                std::cout << "progress : vertex number " << count << " and radius " << r << " out of " << resolution << std::endl; 
                
                dt -= 2 * dummyChord/resolution;
                count += 1;

            }
            r += 1;
        }
        blockMeshGen::generateStl(resolution);
        blockMeshGen::clear();
        //t += 2;
        //}
    }
    
    // blockMeshGen::generateVertices();  
    // blockMeshGen::generateEdges();    
    // blockMeshGen::generateBlocks();
    // blockMeshGen::generateBoundaries(); 

    std::cout << "mesh generation is done\n";
}

/*
void FullNormalEquilibrium()
{

}
*/

void clear()
{
    for(int i = 0; i < 11; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            alpha[i][j].clear();
            beta[i][j].clear();
            Mach[i][j].clear();
            lossCoefficient[i][j].clear();
            liftCoefficient[i][j].clear();
            pressureLoss[i][j].clear();
        }
    }
}

void getAerofoilTD(int i, int j)
{

    char filename[] = "output/stl/blade.stl"; //for blockMesh : blockMeshDict
    blockMeshGen::init(filename);

    double r, shape, extraR;
    r = 1.0;
    shape = 1.0;
    
    fileOut.close();
    // fileOut2.close();
    fileOut.open("output/aerofoil/0.dat" , std::ofstream::trunc);
    // fileOut2.open("aerofoil/circle.dat" , std::ofstream::trunc);
    
    //defining the complex equation
    std::complex<double> z, seta, thetaC, displac, F, aerofoil, dummyF, dummyAerofoil, dummyF1[2], dummyAero1[2];
    double potentialC, streamC, Gamma;

    double multiplier = 1.0;
    double unmultiplier = 1.0;

    //enter the displacement of the circle

    double psiA, psiC, disX, backFat;
    double tempW;

    double U = 1.0;
    psiA = 0.02;
    psiC = 0.012;
    backFat = 1.02;

    double theta;
    double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;
    double radius;
    double dummyScale = chord[i][j];

    double discharge1, discharge2, dischargeX1, dischargeX2;
    int neg = 1;
    int flip = 1;
    
    double threshold = 0.2;
    if(resolution > 50)
    {
    threshold = 0.01 / resolution;
    }

    for(int R = 0; R <= resolution; R++)
    {

    if(j == 0)
    {
        theta = beta[i][0][R] - beta[i][1][R];
        flip = -1;
    }
    if(j == 1)
    {
        theta = alpha[i][1][R] - alpha[i+1][0][R];
        flip = 1;
    }

    radius = hub_radi[i][j] + dr * R;

    disX = psiA - psiC * cos( 0 );
    displac = std::complex( disX , flip * getDisplacementY( liftCoefficient[i][j][R], theta, disX, radius ) );
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );

    dummyF1[0] = std::complex( extraR * cos( 0 )/backFat + displac.real(), extraR * sin(0 )/backFat + displac.imag());
    dummyAero1[0] = joukowskyTransform(dummyF1[0], shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);
    
    disX = psiA - psiC * cos( -PI);
    displac = std::complex( disX , flip * getDisplacementY( liftCoefficient[i][j][R], theta, disX, radius ) );
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );

    dummyF1[1] = std::complex( extraR * cos( -PI )/backFat + displac.real(), extraR * sin(-PI )/backFat + displac.imag());
    dummyAero1[1] = joukowskyTransform(dummyF1[1], shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);
    
    tempW = dummyAero1[0].real() - dummyAero1[1].real();


    while(true)
    {
        for(int l = 0; l < 2;l++)
        {
        
        disX = psiA - psiC * cos( 0 * PI - neg * 1 * PI / resolution );
        displac = std::complex( disX , flip * getDisplacementY( liftCoefficient[i][j][R], theta, disX, radius ) );
        extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );

        dummyF1[l] = std::complex( extraR * cos( 0 * PI - neg * 1 * PI / resolution)/backFat + displac.real(), extraR * sin(0 * PI - 1 * neg * PI / resolution)/backFat + displac.imag());
        dummyAero1[l] = joukowskyTransform(dummyF1[l], shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);
        neg *= -1;
        }

        discharge1 = 0.5 * ( dummyAero1[0].imag() + dummyAero1[1].imag() );
        dischargeX1 = 0.5 * ( dummyAero1[0].real() + dummyAero1[1].real() );
        neg = 1;
        for(int l = 0; l < 2;l++)
        {
        
        disX = psiA - psiC * cos( 0 * PI - neg * 4 * PI / resolution );
        displac = std::complex( disX , flip *  getDisplacementY( liftCoefficient[i][j][R], theta, disX, radius ) );
        extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );

        dummyF1[l] = std::complex( extraR * cos( 0 * PI - neg * 4 * PI / resolution)/backFat + displac.real(), extraR * sin(0 * PI - 4 * neg * PI / resolution)/backFat + displac.imag());
        dummyAero1[l] = joukowskyTransform(dummyF1[l], shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);
        neg *= -1;
        }
        
        discharge2 = 0.5 * ( dummyAero1[0].imag() + dummyAero1[1].imag() );
        dischargeX2 = 0.5 * ( dummyAero1[0].real() + dummyAero1[1].real() );
        
        dischargeAngle[i][j][R] = atan2(  discharge1 - discharge2 ,  dischargeX1 - dischargeX2 ) * RadToDegree;
        
        if((float)beta[i][1][R] == (float)dischargeAngle[i][j][R] or (float)alpha[i+1][0][R] == (float)dischargeAngle[i][j][R] )
        {
            std::cout << "discharge angle matched for stage: " << i << "\n";

            for(double gh = -PI * 1.01 ; gh <= PI * 1.01 ;) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                //fileOut << aerofoil.real() << " " << aerofoil.imag() << std::endl;
                
                blockMeshGen::collectVertices( aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0, radius - hub_radi[i][j] );

                gh += 2.0 * ( PI ) / (double)resolution / multiplier;
            }

            // for(double gh = -PI ; gh <= 0.0;)
            // {
            //     gh += unmultiplier * 2.0 * PI / resolution;

            //     dummyF = std::complex( extraR * cos( gh * 0.9  - 0.2 * (double)resolution) / backFat + displac.real(), extraR * sin(gh * 0.9 - 0.2 * (double)resolution ) / backFat + displac.imag());
            //     aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
            //     blockMeshGen::collectBoundary( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 +  getBladeDist(i, j, radius), radius );
            // }

            // for(double gh = -PI ; gh <= 0.0;)
            // {
            //     gh += unmultiplier * 2.0 * PI / resolution;
                
            //     dummyF = std::complex( extraR * cos( (-PI - gh) * 0.9  - 0.2 * (double)resolution ) / backFat + displac.real(), extraR * sin((-PI - gh) * 0.9  - 0.2 * (double)resolution ) / backFat + displac.imag());
            //     aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
            //     blockMeshGen::collectBoundary( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius );
            // }

            if(R == 0 or R == resolution - 0)
            {
                std::complex<double> dummyG;
                dummyF = std::complex( extraR * cos( -PI ) / backFat + displac.real(), extraR * sin(-PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );

                dummyF = std::complex( extraR * cos( 0.0 ) / backFat + displac.real(), extraR * sin(0.0) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
            
                
                dummyF = std::complex( extraR * cos( -PI ) / backFat + displac.real(), extraR * sin(-PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
            


                dummyF = std::complex( extraR * cos( -0.5 * PI ) / backFat + displac.real(), extraR * sin(-0.5 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );


            }

            // std::cin.get();
            // std::cout << "received\n";

            break;

        }

        if( j == 0)
        {
        rotateAngle[i][j][R] = beta[i][1][R] - dischargeAngle[i][j][R];
        incidenceAngle[i][j][R] = beta[i][0][R] - rotateAngle[i][j][R];
        }
        if( j == 1)
        {
        rotateAngle[i][j][R] = alpha[i+1][0][R] - dischargeAngle[i][j][R];
        incidenceAngle[i][j][R] = alpha[i][1][R] - rotateAngle[i][j][R];
        }
        //std::cout << dischargeAngle[i][j][R] << " " << beta[i][0][R] << std::endl;
        std::cout << "discharge angle error for stage: " << i << " " << R << "\n";

    }

    }
    blockMeshGen::generateStl(multiplier*resolution);
    // blockMeshGen::generateBoundary(resolution/unmultiplier);
    // blockMeshGen::generateInlet(resolution/unmultiplier);
    // blockMeshGen::generateOutlet(resolution/unmultiplier);
    // blockMeshGen::generateBot(resolution/unmultiplier);
    // blockMeshGen::generateTop(resolution/unmultiplier);

    // blockMeshGen::getSeparationVector(getBladeDist(i, j, mean_radi[i][j]));

    blockMeshGen::generateVertices();
    blockMeshGen::generateEdges();
    blockMeshGen::generateBoundaries();
    blockMeshGen::generateBlocks();

    //blockMeshGen::generateBoundaryFile();
    blockMeshGen::generateSnappy();
    blockMeshGen::generateSurfaceFeature();
    blockMeshGen::generateCreatePatch();
    blockMeshGen::generateControlDict();

    blockMeshGen::generateObj(resolution);
    
    
}

void optimiseFlow(int j)    
{
    double theta, yValue, gradient;
    double prevCl = 0.1;
    double dir1, dir2;
    double prevValue;
    int lowestRPM;   

    sVec<int> suitableRPM;
    sVec<double> alphaCombination;

    std::random_device rd1;

    std::uniform_int_distribution<> distr1(50, 400); //omega1
    std::uniform_int_distribution<> distr2(200, 450); //omega2
    std::uniform_int_distribution<> distr3(-60, 60); //alpha1 angles
    std::uniform_int_distribution<> distr4(105, 120); //pressure ratios

    // for(int lda = 0; lda < 100; lda++)
    // {

    // //monte-carlo to find a random combination
    // findCombinationOmega(j,suitableRPM, distr1, distr2, distr3, distr4);

    // //algorithm for finding the lowest rpm
    // //rerolling another random combination
    // omega1 = distr1(rd1) * 20;
    // omega2 = distr2(rd1) * 20;
    // int tempA = distr3(rd1) * 2;
    // double tempB = (double)distr4(rd1) / 100.0; 
    // for(int m = 0; m < lowSize; m++)
    // {
    //     //PR[m] = tempB;
    //     meanAlpha[m][0] = distr3(rd1);

    // }
    // tempA = distr3(rd1) * 2;
    // for(int m = lowSize; m < totalSize; m++)
    // {
    //     //PR[m] = pow( 14 / pow( tempB, 3 ) , 1.0 / 8.0 );
    //     meanAlpha[m][0] = distr3(rd1);
    // }
    // }

    omega1 = findCombinationAlpha(j, alphaCombination, distr2, distr3, 1000);
    for(int b = 0; b < totalSize; b++)
    {
        meanAlpha[b][0] = alphaCombination[b];
    }

    storeOmegaData(alphaCombination);

    //lowestRPM = std::min_element(suitableRPM.begin(), suitableRPM.end())[0];
    
    //storeRandomData();
    
    for(int l = 0; l < totalSize; l++)
    {

        init();
        getFlowPaths(l);

        for(int r = 0; r <= resolution; r++)
        {
            if(j == 0)
            {
            theta = cos(beta[l][0][r]/RadToDegree) / cos(beta[l][1][r]/RadToDegree); 
            }
            if(j == 1)  
            {
            theta = cos(alpha[l+1][0][r]/RadToDegree) / cos(alpha[l][1][r]/RadToDegree); 
            }
            lossCoefficient[l][0][r] = 0.014 * solidity[l][0] / cos( beta[l][1][r] / RadToDegree );
            lossCoefficient[l][1][r] = 0.014 * solidity[l][1] / cos( alpha[l+1][0][r] / RadToDegree );
            
            Mach[l][0][r] = VX[0] / ( cos( beta[l][1][r] / RadToDegree ) * pow( Temperature[l][1] * 287 * gamma , 0.5 ) );
            Mach[l][1][r] = VX[0] / ( cos( alpha[l+1][0][r] / RadToDegree ) * pow( Temperature[l][2] * 287 * gamma , 0.5 ) );

            pressureLoss[l][0][r] = 0.000005 * rho[l][1] * lossCoefficient[l][0][r] * pow( VX[0] / cos( beta[l][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[l][0][r] , 2 ) , gamma / ( gamma - 1 ) );
            pressureLoss[l][1][r] = 0.000005 * rho[l][2] * lossCoefficient[l][1][r] * pow( VX[0] / cos( alpha[l+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[l][1][r] , 2 ) , gamma / ( gamma - 1 ) );

            liftCoefficient[l][0][r] = 2.0 / solidity[l][0] * ( tan( beta[l][0][r] / RadToDegree ) - tan( beta[l][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[l][0][r] / RadToDegree ) + tan( beta[l][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[l][0][r] * sin( 0.5 * ( tan( beta[l][0][r] / RadToDegree ) + tan( beta[l][1][r] / RadToDegree ) ) ) / ( rho[l][1] * pow( 0.5 * ( VX[0] / cos( beta[l][0][r] / RadToDegree ) + VX[0] / cos( beta[l][1][r] / RadToDegree ) ) , 2 ) * solidity[l][0] );
            liftCoefficient[l][1][r] = 2.0 / solidity[l][1] * ( tan( alpha[l][1][r] / RadToDegree ) - tan( alpha[l+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[l][1][r] / RadToDegree ) + tan( alpha[l+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[l][1][r] * sin( 0.5 * ( tan( alpha[l][1][r] / RadToDegree ) + tan( alpha[l+1][0][r] / RadToDegree ) ) ) / ( rho[l][2] * pow( 0.5 * ( VX[1] / cos( alpha[l][1][r] / RadToDegree ) + VX[1] / cos( alpha[l+1][0][r] / RadToDegree ) ) , 2 ) * solidity[l][1] );

            batchAnalysis[l] << theta << " " <<  liftCoefficient[l][j][r] << "\n";
        }
        std::cout<< solidity[l][j] << std::endl;
    }
    std::cout << omega1 << " " << omega2 << std::endl;

}

int getTotalSize()
{
    return totalSize;
}

private:            

template<typename T, typename F>
void insert(dVec<F> vec, T in)
{
    vec.insert(vec.end(), in);
}

void drawShape()
{
    //draw the shape of the compressor, sideways
    int counter = 0;
    for(int i = 0; i < 11; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            compShape << counter << " " << tip_radi[i][j] << "\n";

            counter += 1;
        }
    }
    counter = 0;
    compShape << "\n";
    for(int i = 0; i < 11; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            compShape << counter << " " << hub_radi[i][j] << "\n";

            counter += 1;
        }
    }
}

void drawFlowPath()
{
    double length = 0.0;
    double previousVal = 0.0;
    for(int i = 0; i < 11; i++)
    {
        for(int r = 0; r <= resolution; r++)
        {
            batchAnalysis[i] << length << " " << previousVal + beta[i][0][r] - beta[i][0][0] << "\n";
            length += 1.0;
        }

        previousVal = previousVal + beta[i][0][resolution] - beta[i][0][0];
        batchAnalysis[i] << " " << "\n";

        for(int r = 0; r <= resolution; r++)
        {
            batchAnalysis[i] << length << " " << previousVal + beta[i][1][r] - beta[i][1][0] << "\n";
            length += 1.0;
        }

        previousVal = previousVal + beta[i][1][resolution] - beta[i][1][0];
        batchAnalysis[i] << " " << "\n";

        for(int r = 0; r <= resolution; r++)
        {
            batchAnalysis[i] << length << " " << previousVal + alpha[i][1][r] - alpha[i][1][0] << "\n";
            length += 1.0;
        }

        previousVal = previousVal + alpha[i][1][resolution] - alpha[i][1][0];
        batchAnalysis[i] << " " << "\n";

        for(int r = 0; r <= resolution; r++)
        {
            batchAnalysis[i] << length << " " << previousVal + alpha[i+1][0][r] - alpha[i+1][0][0] << "\n";
            length += 1.0;
        }

        previousVal = previousVal + alpha[i+1][0][resolution] - alpha[i+1][0][0];
        batchAnalysis[i] << " " << "\n";
    }
}

void printOut( std::ofstream &out, double value1, double value2)
{

    out << value1 << " " << value2 << "\n";

}

bool getDeHallerNumber(double alpha1, double alpha2, double beta1, double beta2)
{
    double HallerRotor = cos(alpha2 / RadToDegree )/cos(alpha1 / RadToDegree );
    double HallerStator = cos(beta1 / RadToDegree )/cos(beta2 / RadToDegree );

    std::cout << HallerRotor << " " << HallerStator << std::endl;
    /* if(HallerRotor >= 0.72 && HallerStator >= 0.72 )
    {
        std::cout << "YES" << std::endl;
    } */
    //anthing larger than 0.72 is good
    if( HallerRotor >= 0.72 && HallerStator >= 0.72)
    {
        return 1;
    }
    else
        return 0;
}

void transformAngle(double angle, double &x, double &y)
{
    double dummyX = x;
    double dummyY = y;
    x = dummyX * cos( angle / RadToDegree ) + dummyY * sin( angle / RadToDegree );
    y = -dummyX * sin( angle / RadToDegree ) + dummyY * cos( angle / RadToDegree );
}
//parabolic camberline
// double func(double x, double y) override
// {
//     return ( x * ( dummyChord - x ) ) / ( y * pow( ( dummyChord - 2 * dummyMaxCamPos ) , 2 ) / ( 4 * pow( dummyMaxCam , 2 ) ) + ( dummyChord - 2 * dummyMaxCamPos ) * x / dummyMaxCam - ( ( pow( dummyChord , 2 )  ) - 4 * dummyMaxCamPos * dummyChord ) / ( 4 * dummyMaxCam ) );
// }

//NACA 6-Series camberline
double func(double x, double y) override
{
    // a = dummyMaxCamPos
    // c = dummyChord
    g = -1 / ( 1 - dummyMaxCamPos ) * ( pow( dummyMaxCamPos, 2 ) * ( 0.5 * log( dummyMaxCamPos ) - 0.25 ) + 0.25 );
    h = 1 / ( 1 - dummyMaxCamPos ) * ( 0.5 * pow( 1 - dummyMaxCamPos , 2 ) * log( 1 - dummyMaxCamPos ) - 0.25 * pow( 1 - dummyMaxCamPos , 2 ) ) + g;
    return ( dummyChord * dummyLiftCoefficient / ( 2 * PI * ( dummyMaxCamPos + 1 ) ) * (  1 / ( 1 - dummyMaxCamPos ) * (  0.5 * pow ( dummyMaxCamPos - x / dummyChord , 2 ) * log( fabs( dummyMaxCamPos - x / dummyChord ) ) - 0.5 * pow( 1 - x / dummyChord , 2 ) * log( 1 - x / dummyChord ) + 0.25 * pow( 1 - x / dummyChord , 2 ) - 0.25 * pow( dummyMaxCamPos - x / dummyChord , 2 ) ) - ( x / dummyChord ) * log( x / dummyChord ) + g - h * x / dummyChord ) );

}

double func_real(double x, double y) override
{
    return ( x * ( dummyChord - x ) ) / ( y * pow( ( dummyChord - 2 * dummyMaxCamPos ) , 2 ) / ( 4 * pow( dummyMaxCam , 2 ) ) + ( dummyChord - 2 * dummyMaxCamPos ) * x / dummyMaxCam - ( ( pow( dummyChord , 2 )  ) - 4 * dummyMaxCamPos * dummyChord ) / ( 4 * dummyMaxCam ) );
}

std::complex<double> joukowskyTransform(std::complex<double> z, double shape, std::complex<double> thetaC)
{
    return exp(thetaC) * ( z + pow(shape,2) / z );
}

double getDisplacementY(double Cl, double alpha, double dX, double radius)
{
    return pow(  pow( sin(  asin( -0.5 * Cl / PI  )  ) - alpha / RadToDegree  , 2 ) * ( radius - pow( dX , 2 ) ) / (  1 - pow(  sin(  asin( -0.5 * Cl / PI ) - alpha / RadToDegree )  , 2  )  ) , 0.5 );
}

double getBladeDist(int i, int j, double radius)
{
    double dist = chord[i][j] / solidity[i][j];

    return 0.5 * dist * ( radius / mean_radi[i][j] );
}

void findCombinationOmega(int j, std::vector<int> &out, std::uniform_int_distribution<> distr1,std::uniform_int_distribution<> distr2,std::uniform_int_distribution<> distr3, std::uniform_int_distribution<> distr4 )
{
    float dir1, dir2;
    double theta;
    std::random_device rd1;
    if(j == 1)
    {
        dir1 = 0.8;
        dir2 = 0.0;
    }
    if(j == 0)
    {
        dir1 = 0.0;
        dir2 = -0.8;
    }

    for(int i = 0; i < totalSize ; i++)
    {

    for(int r = 0; r <= resolution; r++)
    {

    while(true)
    {               
    if(j == 0)
    {

    theta = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    }
    if(j == 1)          
    {
    theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    }
    
    lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i][j] / cos( alpha[i+1][0][r] / RadToDegree );
    
    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

    if(theta >= 0.73 && liftCoefficient[i][j][r] <= dir1 && liftCoefficient[i][j][r] >= dir2 && solidity[i][0] <= 5.0)
    {
        std::cout << "stage " << i << " passes\n";      

        break;
    }
    else                                        
    {           
        omega1 = distr1(rd1) * 20;
        omega2 = distr2(rd1) * 20;
        int tempA = distr3(rd1) * 2;
        double tempB = (double)distr4(rd1) / 100.0;
        for(int m = 0; m < lowSize; m++)
        {
            //PR[m] = tempB;
            meanAlpha[m][0] = distr3(rd1);
        }
        tempA = distr3(rd1) * 2;
        for(int m = lowSize; m < totalSize; m++)
        {
            //PR[m] = pow( 14 / pow( tempB, 3 ) , 1.0 / 8.0 );
            meanAlpha[m][0] = distr3(rd1);
        }
        
        std::cout << "trying the following combination: " << omega1 << " " << omega2 << " ";

        init();
        for(int m = 0; m < totalSize; m++)
        {
            std::cout << meanAlpha[m][0] << " ";
            getFlowPaths(m);
        }

        std::cout << "\n";

        r = 0;
        i = 0;
    }

    }

    }
    }

    out.insert(out.end(), omega1);
}

int findCombinationAlpha(int j, sVec<double> &out, std::uniform_int_distribution<> distr2, std::uniform_int_distribution<> distr3, int sampleSize)
{
    float dir1, dir2;
    double theta, smallestGradient;
    std::random_device rd1;
    dVec<int> suitableAlphas;
    dVec<int> successfulAlphas;
    sVec<double> accumulatedSmallestGradient;
    sVec<double> localSmallestGradient;
    sVec<double> stageSmallestGradient;
    dVec<double> tempGradient;
    sVec<double> suitableOmega1;
    int tempOmega1;
    
    tempGradient.reserve(totalSize);
    suitableAlphas.reserve(sampleSize);

    bool firstAttempt = true;
    if(j == 1)
    {
        dir1 = 0.8;
        dir2 = 0.0;
    }
    if(j == 0)
    {
        dir1 = 0.0;
        dir2 = -0.8;
    }
    
    for(int nklf = 0; nklf < sampleSize; nklf++) //sample size
    {

    for(int i = 0; i < totalSize ; i++)
    {

    for(int r = 0; r <= resolution; r++)    
    {
    
    
    while(true)
    {               
    if(j == 0)
    {

    theta = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    }
    if(j == 1)          
    {
    theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    }
    
    lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i][j] / cos( alpha[i+1][0][r] / RadToDegree );
    
    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

    if(theta >= 0.73 && liftCoefficient[i][j][r] <= dir1 && liftCoefficient[i][j][r] >= dir2 && solidity[i][0] <= 5.0 && firstAttempt == false)
    {
        if(r != 0)
        {
        tempGradient[i].insert(tempGradient[i].end(), liftCoefficient[i][j][r] - liftCoefficient[i][j][r-1]);
        }
        if(r == resolution && i == totalSize - 1)
        {
            for( int g = 0; g < totalSize; g++)
            {
                stageSmallestGradient.insert(stageSmallestGradient.end(), std::max_element(tempGradient[g].begin(), tempGradient[g].end())[0]);
                suitableAlphas[nklf].insert(suitableAlphas[nklf].end(), meanAlpha[g][0]);
                suitableOmega1.insert(suitableOmega1.end(), tempOmega1);
                //std::cout << stageSmallestGradient[j] << std::endl;
            }
            localSmallestGradient.insert(localSmallestGradient.end(), std::max_element(stageSmallestGradient.begin(), stageSmallestGradient.end())[0]);
            
            std::cout << "combination no " << nklf << " passes\n";      

            stageSmallestGradient.clear();
            tempGradient.clear();

            firstAttempt = true;
        }
        break;
    }
    else                                        
    {          
        firstAttempt = false; 
        tempGradient[i].clear();     
        suitableAlphas[nklf].clear();

        //int tempA = distr3(rd1) * 2;        
        for(int m = 0; m < lowSize; m++)
        {
            //PR[m] = tempB;
            omega1 = distr2(rd1) * 15;
            tempOmega1 = omega1;
            meanAlpha[m][0] = distr3(rd1);
        }
        //tempA = distr3(rd1) * 2;
        for(int m = lowSize; m < totalSize; m++)
        {
            //PR[m] = pow( 14 / pow( tempB, 3 ) , 1.0 / 8.0 );
            omega1 = distr2(rd1) * 15;
            meanAlpha[m][0] = distr3(rd1);
        }
        
        std::cout << "trying the following combination: " << omega1 << " //" << omega2 << " ";

        init();
        for(int m = 0; m < totalSize; m++)             
        {
            std::cout << meanAlpha[m][0] << " ";
            getFlowPaths(m);
        }

        std::cout << "\n";

        r = 0;
        i = 0;
    }

    }

    }
    }

    //accumulatedSmallestGradient.insert(accumulatedSmallestGradient.end(), std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0]);

    }
    // for(int l = 0; l < sampleSize; l++)
    // {
    //     std::cout << localSmallestGradient[l] << std::endl;
    // }
    smallestGradient = std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0];
    auto index = std::find(localSmallestGradient.begin(),localSmallestGradient.end(), smallestGradient);
    size_t dist = std::distance(localSmallestGradient.begin(), index);

    std::cout << dist << std::endl;

    for(int c = 0; c < totalSize; c++)
    {
    out.insert(out.end(), suitableAlphas[(int)dist][c]);
    }
    return tempOmega1;
}


void storeRandomData()
{
    std::ofstream lk("output/randomGenResult/angles.dat", std::ios::app);
    std::ofstream lo("output/randomGenResult/omegas.dat", std::ios::app);


    for(int l = 0; l < totalSize; l++)
    {
        lk << l << " " << meanAlpha[l][0] << std::endl;     
    }

    lk << "\n";

    lo << omega1 << " " << omega2 << "\n";

    lk.close();
    lo.close();
};

void storeOmegaData(sVec<double> omegas)
{
    std::ofstream lk("output/randomGenResult/angles.dat", std::ios::app);
    std::ofstream lo("output/randomGenResult/omegas.dat", std::ios::app);

    lo << omega1 << "\n";

    lo.close();

    for(int l = 0; l < totalSize; l++)
    {
        lk << l << " " << omegas[l] << std::endl;     
    }

    lk << "\n";

    lk.close();
};

};

int main()
{
    dVec<double> rTip = {
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },           
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    { 0.08, 0.08, 0.08 },
    };
    sVec<double> delta_P = { 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}; //{ 1.075, 1.1, 1.15, 1.2, 1.25, 1.3, 1.3, 1.325, 1.325, 1.325, 1.3 };
    sVec<double> init_alpha1 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }; //{ -50, -48, -46, 0, 30, 42, 48, 58, 58, 62, 68, 70 };
    sVec<double> degOfReaction = { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,0.6, 0.6};
    dVec<double> chordLengths = {
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },          
        { 0.015 , 0.015 },
        { 0.015 , 0.015 },  
        { 0.015 , 0.015 },          
    };                          
    double rHub = 0.04;
    double omega_1 = 2250;//5390; //3250 //from algorithm = 3475
    double omega_2 = 5600;//6650; //5600 //from algorithm = 7275
    //resolution only works with even numbers, idk why
    double resol = 8; 
    double v1 = 69.4377418988;                                              
    double v2 =  v1;
    double Temp = 300;
    double rho_ = 1.204;
    double Pres = 103.15;

    //system("./lib/clear.sh");

    Blade test(rTip, rHub, omega_1, omega_2 , resol, v1, v2, delta_P, Temp, Pres, init_alpha1, degOfReaction, chordLengths, 7, 0);

    test.init();
    
    int f = 0;
    for(int i = 0; i < test.getTotalSize(); i++)
    {
    test.getFlowPaths(i);
    }
    //test.analysisCamber(f);              
    test.optimiseFlow(f);
    //test.getLeadingTrailingAngles(f);
    // for(int i = 0; i < 11; i++)                  
    // {                        
    // test.getBladeAngles(i,f);
    // }
    //test.getAerofoilTD(3,f);
    //test.generateBlade(0,f);
    //test.getCamberline(3,1,50);
    //test.clear();


    return 0;
}

//plotting with gnuplot:
//plot '<datafile.dat>' with linespoints linetype 0 linewidth 2
//plot 'out.dat' with linespoints linetype 0 linewidth 2, 'out2.dat' with linespoints linetype 0 linewidth 2, 'out3.dat' with linespoints linetype 0 linewidth 2, 'out4.dat' with linespoints linetype 0 linewidth 2
//compile : g++ main.cpp ODE_solver/solver1.o -o EXE
//gnuplot 1:1 ratio : set size ratio -1
//gnuplot camberline :  set size ratio -1
//                      plot'camberline.dat' with linespoints linetype 0 linewidth 2
//g++ ODE_solver/solver1.cpp blockMeshGenerator/blockMeshGenerator.cpp main.cpp -o EXE && ./EXE


//{0.09, 0.09, 0.089 },
// { 0.089, 0.088, 0.087 },
// { 0.087, 0.086, 0.084 },
// { 0.084, 0.02, 0.081 },
// { 0.081, 0.08, 0.08 },
// { 0.08, 0.08, 0.08 },
// { 0.08, 0.08, 0.08 },
// { 0.08, 0.08, 0.08 },           
// { 0.08, 0.08, 0.08 },
// { 0.08, 0.08, 0.08 },
// { 0.08, 0.08, 0.08 },