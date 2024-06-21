#include <filesystem>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <ratio>
#include <thread>
#include <vector>
#include "ODE_solver/solver1.h"
#include "blockMeshGenerator/blockMeshGenerator.h"

#define e 2.718281828459045
#define PI 3.14159265
#define RadToDegree 57.29577951
#define Cp 1004.5
#define gamma 1.4

class Blade : public solver1, public blockMeshGen
{

std::ofstream fileOut;
std::ofstream fileOut2;
std::ofstream fileOut3;
std::ofstream fileOut4;

std::ofstream batchAnalysis[11];

double VX[3];
double tip_radi[11][3];
double hub_radi[11][3];
double mean_radi[11][3];
double omega1;
double omega2;
double v[11][3];
//resolution is 100 here 
double alpha[12][3][200];
double beta[12][2][200];
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
//dummy variables
double dummyMaxCam, dummyMaxCamPos, dummyChord;
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

double solidity[11];

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

        }

        //opening a file for each batchAnalysis ofstream
        tempName = "batchData/" + std::to_string(i) + ".dat";
        batchAnalysis[i].open( tempName);
    }

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
    double work = Cp * ( TemperatureStag[0][2] - TemperatureStag[0][0] );
    TemperatureStag[0][1] = TemperatureStag[0][0] + ( R[0] * work / Cp );
    //Pressure 2
    Pressure[0][1] = Pressure[0][0] * pow( ( TemperatureStag[0][1] / TemperatureStag[0][0] ) , gamma / ( gamma - 1 ) );
    //Stagnation pressure 2
    PressureStag[0][1] = PressureStag[0][0] * pow( ( TemperatureStag[0][1] / TemperatureStag[0][0] ) , gamma / ( gamma - 1 ) );
    //psi
    psi[0] = work / pow( omega1 * mean_radi[0][0] , 2 );
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
    double vu2_r = VX[0] * tan( meanAlpha[0][1] / RadToDegree );
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
    double work = Cp * ( TemperatureStag[i][2] - TemperatureStag[i][0] );
    TemperatureStag[i][1] = TemperatureStag[i][0] + ( R[i] * work / Cp );
    //Pressure 2
    Pressure[i][1] = Pressure[i][0] * pow( ( TemperatureStag[i][1] / TemperatureStag[i][0] ) , gamma / ( gamma - 1 ) );
    //Stagnation pressure 2
    PressureStag[i][1] = PressureStag[i][0] * pow( ( TemperatureStag[i][1] / TemperatureStag[i][0] ) , gamma / ( gamma - 1 ) );
    //psi
    if( i < 3 )
    {
    psi[i] = work / pow( omega1 * mean_radi[i][0] , 2 );
    }
    if( i >= 3 )
    {
    psi[i] = work / pow( omega2 * mean_radi[i][0] , 2 );
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
    double vu2_r = VX[0] * tan( meanAlpha[i][1] / RadToDegree );
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

        tempVu1 = a[i] - b[i]/y;
        tempVu2 = a[i] + b[i]/y;
        
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

        //printOut(fileOut, r, alpha[i][0][r]); 
        //printOut(fileOut2, r, alpha[i][1][r]);
        //printOut(fileOut3, r, beta[i][0][r]);
        //printOut(fileOut4, r, beta[i][1][r]);       
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
    double incidence = 0;
    //chosen stagger angle
    double stagger = 35;
    double Ksh = 1.0;
    double tb = 0.1 * chord[i][0];
    //the distance to maximum thickness, might need to change it later
    double distToMax = 0.5;
    bool stop = 0;
    int correctionVar = 1;
    //everything is positive, need to readjust them again later
    for(int r = 0; r <= resolution; r++)
    {
        stop = 0;
        radius1 = hub_radi[i][0] + r * dr1;
        radius2 = hub_radi[i][1] + r * dr2;
        k1 = beta[i][0][r];
        k2 = beta[i][1][r];
        dummyK1 = k1;
        dummyK2 = k2;

        //need a loop for the following, use the stop boolean and if statement commented below
        if(j == 0)
        {
        theta = beta[i][0][r] - beta[i][1][r]; 
        }
        if(j == 1)
        {
        theta = alpha[i][1][r] - alpha[i+1][0][r]; 
        }
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

        delta0Point = 0.01 * solidity[i] * beta[i][0][r] + ( 0.74 * pow( solidity[i] , 1.9 ) + 3 * solidity[i] ) * pow( fabs( beta[i][0][r] ) / 90 , 1.67 + 1.09 * solidity[i] );
        tempX = beta[i][0][r] / 100;
        m10 = 0.17 - 0.0333 * tempX + 0.333 * pow( tempX , 2 );
        tempB = 0.9625 - 0.17 * tempX - 0.85 * pow( tempX , 3 );
        m = m10 / pow( solidity[i] , tempB );
        Ktdelta = 6.25 * ( tb / chord[i][0] ) + 37.5 * pow( tb / chord[i][0] , 2 );
        deltaPoint = Ksh * Ktdelta * delta0Point + m * theta;

        dIncidence = incidence - incidencePoint;
        //TODO fix incidencePoint
        k1 = beta[i][0][r] - incidencePoint - dIncidence;
        k2 = beta[i][1][r] - deltaPoint;
        
        std::cout << k1 << " " << alphaPoint << std::endl;
        
        /* //convergence test
        if ( k1 / 1.05 <= dummyK1 <= k1 * 1.05 && k2 / 1.05 <= dummyK2 <= k2 * 1.05 )
        {
            stop = 1;
        } */

    }
}

//j is either one(rotor) or two(stator)
void getCamberline(int i, int j, int r)
{
    double dummyR = 0.0;
    //distance to the maximum camber, dimesionless
    maxCamPos[i][j] = 0.42 * chord[i][j];
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

// j = 0 for rotor and j = 1 for stator
void analysisCamber(int j)
{
    for(int i = 0; i < 11; i++)
    {
    //maxCamPos[i][j] = 0.4 * chord[i][j];
    double theta, yValue;

    //double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;

    for(int r = 0; r <= resolution; r++)
    {
    //yValue = ( hub_radi[i][j] + dr * r ) / mean_radi[i][j];
    // if(j == 0)
    // {
    // theta = cos (beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    // }
    // if(j == 1)
    // {
    // theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    // }

    if(j == 0)
    {
    theta = beta[i][0][r] - beta[i][1][r]; 
    }
    if(j == 1)
    {
    theta = alpha[i][1][r] - alpha[i+1][0][r]; 
    }

    // maxCam[i][j] = chord[i][j] / ( 4 * tan( theta / RadToDegree ) ) * ( sqrt( fabs( 1 + pow( 4 * ( tan( theta / RadToDegree ) ) , 2 ) * ( maxCamPos[i][j] / chord[i][j] - pow( maxCamPos[i][j] / chord[i][j] - 3.0 / 16.0 , 2 ) ) ) ) - 1 );

    // if(j == 1)
    // {
    // batchAnalysis[i] << r << " " << alpha[i][j][r] << "\n";
    // }
    // if(j == 0)
    // {
    // batchAnalysis[i] << r << " " << beta[i][j][r] << "\n";
    // }
    batchAnalysis[i] << r << " " <<  theta << "\n";
    
    }
    }
    //solver1::rungeKuttam(dummyChord, 0.0, 0.01 / resolution, 0.0, dummyChord, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);

}

void generateBlade(int i, int j)
{
    double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;
    double radius;

    char filename[] = "blade.stl"; //for blockMesh : blockMeshDict
    blockMeshGen::init(filename);

    for(int r = 0; r <= resolution; r++)
    {
        radius = dr * r + hub_radi[i][j];

        double dummyR = 0.0;
        //distance to the maximum camber, dimesionless
        maxCamPos[i][j] = 0.42 * chord[i][j];
        double theta;

        if(j == 0)
        {
        theta = beta[i][0][r] - beta[i][1][r]; 
        }
        if(j == 1)
        {
        theta = alpha[i][1][r] - alpha[i+1][0][r]; 
        }

        maxCam[i][j] = chord[i][j] / ( 4 * tan( theta / RadToDegree ) ) * ( sqrt( fabs( 1 + pow( 4 * ( tan( theta / RadToDegree ) ) , 2 ) * ( maxCamPos[i][j] / chord[i][j] - pow( maxCamPos[i][j] / chord[i][j] - 3.0 / 16.0 , 2 ) ) ) ) - 1 );

        //initialising dummy variables
        for(int t = -1; t < 2;)
        {
        dummyMaxCam = maxCam[i][j] + t * 0.05 * chord[i][j];
        dummyMaxCamPos = maxCamPos[i][j];
        dummyChord = chord[i][j]; 


        //std::cout << theta << " " << dummyMaxCam << " " << dummyMaxCamPos << " " << dummyChord << std::endl;

        int count = 0;
        for(double dt = 0; dt <= dummyChord;)
        {
            dummyR = func( dt, dummyR );
            dt += dummyChord/resolution;

            blockMeshGen::collectVertices(dt, dummyR, radius);

            std::cout << "progress : vertex number " << count << " and radius " << r << " out of " << resolution << std::endl; 
            count += 1;
        }
        t += 2;
        }
    }

    blockMeshGen::generateStl(resolution);
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

private:

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

double func(double x, double y) override
{
    return ( x * ( dummyChord - x ) ) / ( y * pow( ( dummyChord - 2 * dummyMaxCamPos ) , 2 ) / ( 4 * pow( dummyMaxCam , 2 ) ) + ( dummyChord - 2 * dummyMaxCamPos ) * x / dummyMaxCam - ( ( pow( dummyChord , 2 )  ) - 4 * dummyMaxCamPos * dummyChord ) / ( 4 * dummyMaxCam ) );
}

double func_real(double x, double y) override
{
    return ( x * ( dummyChord - x ) ) / ( y * pow( ( dummyChord - 2 * dummyMaxCamPos ) , 2 ) / ( 4 * pow( dummyMaxCam , 2 ) ) + ( dummyChord - 2 * dummyMaxCamPos ) * x / dummyMaxCam - ( ( pow( dummyChord , 2 )  ) - 4 * dummyMaxCamPos * dummyChord ) / ( 4 * dummyMaxCam ) );
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
    double delta_P[11] = { 1.125, 1.15, 1.15, 1.225, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25 };
    double init_alpha1[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 10, 10 };
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
    double omega_1 = 6750;
    double omega_2 = 9850;
    double resol = 50;
    double v1 = 69.4377418988;
    double v2 =  v1;
    double Temp = 300;
    double rho_ = 1.204;
    double Pres = 103.15;

    Blade test(rTip, rHub, omega_1, omega_2 , resol, v1, v2, delta_P, Temp, Pres, init_alpha1, degOfReaction, chordLengths);

    test.init();

    //angles are not available yet, need to initialise them for every r and stage
    for(int i = 0; i < 11; i++)
    {
    test.getFlowPaths(i);
    }
    //test.getBladeAngles(0,0);
    test.generateBlade(10,1);
    //test.getCamberline(3,1,50);
    //test.analysisCamber(1);

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