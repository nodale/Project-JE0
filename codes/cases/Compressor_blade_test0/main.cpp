#include <filesystem>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <ratio>
#include <chrono>
#include <thread>
#include <vector>
#include "ODE_solver/solver1.h"
#include "blockMeshGenerator/blockMeshGenerator.h"
#include <complex>

#define e 2.718281828459045
#define PI 3.14159265
#define RadToDegree 57.29577951
#define Cp 1004.5
#define gamma 1.4

//TODO
//finish conformal transform, psi, ehta, and y are variables of x

class Blade : public solver1, public blockMeshGen
{

std::ofstream fileOut;
std::ofstream fileOut2;
std::ofstream fileOut3;
std::ofstream fileOut4;

std::ofstream compShape;

std::ofstream batchAnalysis[11];

double VX[3];
double tip_radi[11][3];
double hub_radi[11][3];
double mean_radi[11][3];
double omega1;
double omega2;
double v[11][3];
double work[11];

//Mach[i][0] is for velocity 2 relative
//Mach[i][1] is for velocity 3 absolute 
std::vector<double> Mach[11][2];

std::vector<double> alpha[12][2];
std::vector<double> beta[12][2];
double meanAlpha[12][2];
double meanBeta[11][2];
double Mach0;
double resolution;
double PR[11];
double R[11];
double Area[11][3];
double chord[11][2];
double maxCam[11][2];
double maxCamPos[11][2];
std::vector<double> liftCoefficient[11][2];
std::vector<double> rotateAngle[11][2];
std::vector<double> incidenceAngle[11][2];
std::vector<double> AoA[11][2];
//dummy variables
double dummyMaxCam, dummyMaxCamPos, dummyChord, dummyLiftCoefficient;
//work-done factor
double WDF[18] = { 0.982, 0.952, 0.929, 0.910, 0.895, 0.882, 0.875, 0.868, 0.863, 0.860, 0.857, 0.855, 0.853, 0.851, 0.850, 0.849, 0.848, 0.847 };

//thermodynamics initial conditions
double T1;
double rho1;
double P1;

double Temperature[11][3];
double TemperatureStag[11][3];
double Pressure[11][3];
double PressureStag[11][3];
double rho[11][3];
double psi[11];
double phi[11];
double a[11];
double b[11];
double Wr[11];
double Ws[11];
double efficiency[11][2];
std::vector<double> lossCoefficient[11][2];
std::vector<double> pressureLoss[11][2];
std::vector<double> dischargeAngle[11][2];
double aValues[5] = { 0.2969, -0.1260, -0.3516, 0.2843, -0.1015 } ;

double solidity[11];


double h, g;

public:
Blade(double (&tipP)[11][3], double hub, double w1, double w2, double res, double vx1, double vx2, double (&deltaP)[11], double Temp, double P, double (&initial_alpha1)[12], double (&Reaction)[11], double (&c)[11][2])
: omega1 { w1 }
, omega2 { w2 }
, resolution { res }
, VX { vx1, vx2, vx1 }
, T1 { Temp }
, P1 { P }
{
    fileOut.open("out.dat");
    fileOut2.open("out2.dat");
    fileOut3.open("out3.dat");
    fileOut4.open("out4.dat");
    compShape.open("shape.dat");
    char camFilename[] = "camberline.dat";
    solver1::initFile(camFilename);

    std::string tempName;
    for(int i = 0; i < 11; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            tip_radi[i][j] = tipP[i][j];

            PR[i] = deltaP[i];

            meanAlpha[i][0] = initial_alpha1[i];

            R[i] = Reaction[i];

            chord[i][0] = c[i][0];
            chord[i][1] = c[i][1];
            alpha[i][j].reserve(resolution);
            beta[i][j].reserve(resolution);
            lossCoefficient[i][j].reserve(resolution);
            liftCoefficient[i][j].reserve(resolution);
            Mach[i][j].reserve(resolution);
            pressureLoss[i][j].reserve(resolution);
            dischargeAngle[i][j].reserve(resolution);
            rotateAngle[i][j].reserve(resolution);
            incidenceAngle[i][j].reserve(resolution);
            AoA[i][j].reserve(resolution);
            //this doesn't work for some reasons
            //efficiency[i][j] = 0.9;
            
        }

        //opening a file for each batchAnalysis ofstream
        tempName = "batchData/" + std::to_string(i) + ".dat";
        batchAnalysis[i].open( tempName);
    }
    
    meanAlpha[11][0] = initial_alpha1[11];

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

    solidity[0] = ( 1.5 * psi[0] ) / ( 1.55 * phi[0] - psi[0] );

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
    for(int i = 1; i < 11; i++)
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
    if( i < 3 )
    {
    phi[i] =  VX[0] / ( omega1 * mean_radi[i][0] );
    }
    if( i >= 3 )
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
    if( i < 3 )
    {
    psi[i] = work[i] / pow( omega1 * mean_radi[i][0] , 2 );
    }
    if( i >= 3 )
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

    solidity[i] = ( 1.5 * psi[i] ) / ( 1.55 * phi[i] - psi[i] );

    double vu1_r = VX[0] * tan( meanAlpha[i][0] / RadToDegree );
    // approach 1
    //double vu2_r = VX[0] * tan( meanAlpha[i][1] / RadToDegree );
    
    //approach 2
    double vu2_r;
    if(i < 3)
    {
        vu2_r = ( work[i] + omega1 * vu1_r * mean_radi[i][0] ) / ( omega1 * mean_radi[i][0] );
    }
    if(i >= 3)
    {
        vu2_r = ( work[i] + omega2 * vu1_r * mean_radi[i][0] ) / ( omega2 * mean_radi[i][0] );
    }

    
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
        if( i < 3)
        {
        tempPhi = VX[0] / ( omega1 * radius1 );
        }
        if(i >= 3)
        {
        tempPhi = VX[0] / ( omega2 * radius1 );  
        }

        alpha[i][0][r] = atan( tempVu1 / VX[0] ) * RadToDegree;
        alpha[i][1][r] = atan( tempVu2 / VX[0] ) * RadToDegree;

        beta[i][0][r] = atan( tan( alpha[i][0][r] / RadToDegree ) - 1 / tempPhi ) * RadToDegree;
        beta[i][1][r] = atan( tan( alpha[i][1][r] / RadToDegree ) - 1 / tempPhi ) * RadToDegree;

        alpha[11][0][r] = meanAlpha[11][0];

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

// j = 0 for rotor and j = 1 for stator, can be used for any analysis
void analysisCamber(int j)
{
    for(int i = 0; i < 11; i++)
    {
    maxCamPos[i][j] = 0.5 * chord[i][j];
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
    lossCoefficient[i][0][r] = 0.014 * solidity[i] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i] / cos( alpha[i+1][0][r] / RadToDegree );
    
    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    //pressure loss is in Kpascal
    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[i] / cos( beta[i][0][r] / RadToDegree ) + VX[i] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[i] / cos( alpha[i][1][r] / RadToDegree ) + VX[i] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i] );

    batchAnalysis[i] << theta << " " <<  liftCoefficient[i][j][r] << "\n";

    //efficiency analysis
    //double dummyEfficiency = 1 - ( lossCoefficient[i][0][r] / pow ( cos( alpha[i+1][0][r] / RadToDegree ), 2 ) + lossCoefficient[i][1][r] / pow ( cos( beta[i][0][r] / RadToDegree ), 2 ) ) * pow( phi[i] , 2 ) / ( 2 * psi[i] );
    //batchAnalysis[i] << r << " " <<  dummyEfficiency << "\n";
    }
    }
    //solver1::rungeKuttam(dummyChord, 0.0, 0.01 / resolution, 0.0, dummyChord, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);

    //drawShape();
    getAerofoilTD();
}

void getLeadingTrailingAngles(int j)
{
    double dt, discharge1, discharge2, dummyDischargeAngle, dummyR;
    for(int i = 0; i < 11; i++)
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
            delta0Point = 0.01 * solidity[i] * beta[i][0][r] + ( 0.74 * pow( solidity[i] , 1.9 ) + 3 * solidity[i] ) * pow( fabs( beta[i][0][r] ) / 90 , 1.67 + 1.09 * solidity[i] );
            }
            if(j == 1)
            {
            delta0Point = 0.01 * solidity[i] * beta[i][0][r] + ( 0.74 * pow( solidity[i] , 1.9 ) + 3 * solidity[i] ) * pow( fabs( alpha[i][1][r] ) / 90 , 1.67 + 1.09 * solidity[i] );
            }

            tempX = beta[i][0][r] / 100;
            m10 = 0.17 - 0.0333 * tempX + 0.333 * pow( tempX , 2 );
            tempB = 0.9625 - 0.17 * tempX - 0.85 * pow( tempX , 3 );
            m = m10 / pow( solidity[i] , tempB );
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

void getAerofoilTD()
{

    double r, shape, psiA;
    r = 1.0;
    shape = 0.02;
    psiA = 0.12;
    fileOut.close();
    fileOut.open("aerofoil/0.dat" , std::ofstream::trunc);
    

    //defining the complex equation
    std::complex<double> z, seta, thetaC;

    //getting the abscissa and ordinates

    for( double thetaAerofoil = 0; thetaAerofoil <= 2 * PI;)
    {
        thetaC = std::complex<double>( 0.0 , 45 / RadToDegree);
        z = std::complex( r * ( exp(psiA) + exp(-psiA) ) * cos(thetaAerofoil), r * ( exp(psiA) - exp(-psiA) ) * sin(thetaAerofoil) );
        seta = joukowskyTransform(z, shape);



        thetaAerofoil += PI / resolution;
        fileOut << seta.real() << " " << seta.imag() << "\n";
    }
    // xAerofoil = 2 * scaleAerofoil * coshl( psiAerofoil ) * cos( thetaAerofoil );
    // yAerofoil = 2 * scaleAerofoil * sinhl( psiAerofoil ) * sin( thetaAerofoil );

    // thetaAerofoil
    // yAerofoil = 2 * scaleAerofoil * sinhl( psiAerofoil ) * sin( thetaAerofoil );

    
}

private:

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

std::complex<double> joukowskyTransform(std::complex<double> z, double shape)
{
    return ( z + pow(shape,2) / z );
}

};

int main()
{
    double rTip[11][3] = {
        {0.09, 0.09, 0.089 },
        { 0.089, 0.088, 0.087 },
        { 0.087, 0.086, 0.081 },
        { 0.081, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
        { 0.08, 0.08, 0.08 },
    };
    double delta_P[11] = { 1.075, 1.1, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.325, 1.325, 1.325 };
    double init_alpha1[12] = { -50, -48, -46, 0, 30, 42, 48, 58, 58, 62, 68, 70 };
    double degOfReaction[11] = { 0.5, 0.5, 0.5, 0.6, 0.6, 0.6,0.6, 0.6, 0.6, 0.6, 0.6 };
    double chordLengths[11][2] = {
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
    double rHub = 0.05;
    double omega_1 = 3250; //3250
    double omega_2 = 5600; //5600
    //resolution only works best with X50; where X is any integer
    double resol = 150; 
    double v1 = 69.4377418988;
    double v2 =  v1;
    double Temp = 300;
    double rho_ = 1.204;
    double Pres = 103.15;

    Blade test(rTip, rHub, omega_1, omega_2 , resol, v1, v2, delta_P, Temp, Pres, init_alpha1, degOfReaction, chordLengths);

    test.init();
    
    //angles are not available yet, need to initialise them for every r and stage
    int f = 0;
    for(int i = 0; i < 11; i++)
    {
    test.getFlowPaths(i);
    }
    test.analysisCamber(f);
    //test.getLeadingTrailingAngles(f);
    // for(int i = 0; i < 11; i++)
    // {
    // test.getBladeAngles(i,f);
    // }
    //test.generateBlade(0,f);
    //test.getCamberline(3,1,50);
    //test.clear();
    
    //TODO check De Haller's number

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