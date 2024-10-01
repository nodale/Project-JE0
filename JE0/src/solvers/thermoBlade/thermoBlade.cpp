#include "thermoBlade.h"

void readTipRadius(dVec<double>& tipRadius)
{
    std::ifstream input("input/tipRadiusConfig.dat");
    std::string temp;
    double tempD;
    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        std::getline(input, temp);
        std::istringstream lineStream(temp);
        for(int j = 0; j < 3; j++)
        {
            lineStream >> tempD;
            tipRadius[i].insert(tipRadius[i].begin() + j, tempD);
        }
    }
    input.close();
}

void readChordLengths(dVec<double>& chordLengths)
{
    std::ifstream input("input/chordConfig.dat");
    std::string temp;
    double tempD;
    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        std::getline(input, temp);
        std::istringstream lineStream(temp);
        for(int j = 0; j < 2; j++)
        {
            lineStream >> tempD;
            chordLengths[i][j] = tempD;
        }
    }
    input.close();
}

void readAngles(dVec<double>& meanAlpha)
{
    std::ifstream input("input/initialAlphaConfig.dat");
    std::string temp;
    for(int i = 0; i <= infoBlade::totalSize; i++)
    {
        std::getline(input, temp);
        meanAlpha[i][0] = std::stod(temp);
    }
    input.close();
}

void readInputFile(sVec<double>& inputMatrix, const char* argc)
{

    std::ifstream input(argc);
    std::string temp;
    int i = 0;
    while(std::getline(input, temp))
    {
        inputMatrix[i] = std::stod(temp);
        i++;
    }

    input.close();
}

void storeInDatabaseRecursive(sqlite3* db)
{
    using namespace infoBlade;
    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        storeInThermoDatabase(db, "Temperature1", Temperature[i][0], i + 1);
        storeInThermoDatabase(db, "Temperature2", Temperature[i][1], i + 1);
        storeInThermoDatabase(db, "Temperature3", Temperature[i][2], i + 1);

        storeInThermoDatabase(db, "TemperatureStag1", TemperatureStag[i][0], i + 1);
        storeInThermoDatabase(db, "TemperatureStag2", TemperatureStag[i][1], i + 1);
        storeInThermoDatabase(db, "TemperatureStag3", TemperatureStag[i][2], i + 1);

        storeInThermoDatabase(db, "Pressure1", Pressure[i][0], i + 1);
        storeInThermoDatabase(db, "Pressure2", Pressure[i][1], i + 1);
        storeInThermoDatabase(db, "Pressure3", Pressure[i][2], i + 1);

        storeInThermoDatabase(db, "PressureStag1", PressureStag[i][0], i + 1);
        storeInThermoDatabase(db, "PressureStag2", PressureStag[i][1], i + 1);
        storeInThermoDatabase(db, "PressureStag3", PressureStag[i][2], i + 1);
        
        storeInThermoDatabase(db, "rho1", rho[i][0], i + 1);
        storeInThermoDatabase(db, "rho2", rho[i][1], i + 1);
        storeInThermoDatabase(db, "rho3", rho[i][2], i + 1);

        storeInThermoDatabase(db, "psi", psi[i], i + 1);

        storeInThermoDatabase(db, "phi", phi[i], i + 1);

        storeInThermoDatabase(db, "a", a[i], i + 1);
        storeInThermoDatabase(db, "b", b[i], i + 1);

        storeInThermoDatabase(db, "Wr",Wr[i], i + 1);
        storeInThermoDatabase(db, "Ws",Ws[i], i + 1);

        storeInThermoDatabase(db, "efficiency1",efficiency[i][0], i + 1);
        storeInThermoDatabase(db, "efficiency2",efficiency[i][1], i + 1);
        //storeInDatabase(db, "efficiency3",efficiency[i][2], i + 1);
    }
}

void resizeAndReserve()
{
    infoBlade::fileOut.open("output/misc/out.dat");    
    infoBlade::fileOut2.open("output/misc/out2.dat");   
    infoBlade::fileOut3.open("output/misc/out3.dat");  
    infoBlade::fileOut4.open("output/misc/out4.dat");
    infoBlade::compShape.open("output/misc/shape.dat");
    char camFilename[] = "output/misc/camberline.dat";

    infoBlade::tip_radi.resize(infoBlade::totalSize);
    infoBlade::hub_radi.resize(infoBlade::totalSize);
    infoBlade::mean_radi.resize(infoBlade::totalSize);
              
    infoBlade::v.resize(infoBlade::totalSize);
    infoBlade::work.resize(infoBlade::totalSize);
    infoBlade::diffusion.resize(infoBlade::totalSize);
    infoBlade::numBlades.resize(infoBlade::totalSize);

    infoBlade::Mach.resize(infoBlade::totalSize);

    infoBlade::alpha.resize(infoBlade::totalSize + 1);
    infoBlade::beta.resize(infoBlade::totalSize + 1);
    infoBlade::meanAlpha.resize(infoBlade::totalSize + 1);
    infoBlade::meanBeta.resize(infoBlade::totalSize);
    infoBlade::PR.resize(infoBlade::totalSize);
    infoBlade::R.resize(infoBlade::totalSize);
    infoBlade::Ry.resize(infoBlade::totalSize);
    infoBlade::Area.resize(infoBlade::totalSize + 1);
    infoBlade::chord.resize(infoBlade::totalSize);
    infoBlade::maxCam.resize(infoBlade::totalSize);       
    infoBlade::maxCamPos.resize(infoBlade::totalSize);
    infoBlade::liftCoefficient.resize(infoBlade::totalSize);
    infoBlade::rotateAngle.resize(infoBlade::totalSize);
    infoBlade::incidenceAngle.resize(infoBlade::totalSize);
    infoBlade::AoA.resize(infoBlade::totalSize);

    infoBlade::Temperature.resize(infoBlade::totalSize);
    infoBlade::TemperatureStag.resize(infoBlade::totalSize);
    infoBlade::Pressure.resize(infoBlade::totalSize);
    infoBlade::PressureStag.resize(infoBlade::totalSize);
    infoBlade::rho.resize(infoBlade::totalSize);
    infoBlade::psi.resize(infoBlade::totalSize);
    infoBlade::phi.resize(infoBlade::totalSize);
    infoBlade::a.resize(infoBlade::totalSize);
    infoBlade::b.resize(infoBlade::totalSize);
    infoBlade::Wr.resize(infoBlade::totalSize);
    infoBlade::Ws.resize(infoBlade::totalSize);
    infoBlade::efficiency.resize(infoBlade::totalSize);          
    infoBlade::lossCoefficient.resize(infoBlade::totalSize);
    infoBlade::pressureLoss.resize(infoBlade::totalSize);
    infoBlade::dischargeAngle.resize(infoBlade::totalSize);

    infoBlade::solidity.resize(infoBlade::totalSize);

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        infoBlade::hub_radi[i].resize(3);
        infoBlade::mean_radi[i].resize(3);
        infoBlade::v[i].resize(3);
        infoBlade::diffusion[i].resize(2);
        infoBlade::numBlades[i].resize(2);
        infoBlade::Ry[i].resize(infoBlade::resolution);

        infoBlade::Mach[i].resize(2);

        infoBlade::alpha[i].resize(2);
        infoBlade::beta[i].resize(2);
        infoBlade::meanAlpha[i].resize(2);
        infoBlade::meanBeta[i].resize(2);
        infoBlade::Area[i].resize(3);
        infoBlade::chord[i].resize(2);
        infoBlade::maxCam[i].resize(2);       
        infoBlade::maxCamPos[i].resize(2);
        infoBlade::liftCoefficient[i].resize(2);
        infoBlade::rotateAngle[i].resize(2);
        infoBlade::incidenceAngle[i].resize(2);
        infoBlade::AoA[i].resize(2);
        
        infoBlade::Temperature[i].resize(3);
        infoBlade::TemperatureStag[i].resize(3);
        infoBlade::Pressure[i].resize(3);
        infoBlade::PressureStag[i].resize(3);
        infoBlade::rho[i].resize(3);
        infoBlade::efficiency[i].resize(2);       
        infoBlade::lossCoefficient[i].resize(2);
        infoBlade::pressureLoss[i].resize(2);
        infoBlade::dischargeAngle[i].resize(2);
        
        infoBlade::solidity[i].resize(2);
    }

    infoBlade::meanAlpha[infoBlade::totalSize].resize(2);
    infoBlade::alpha[infoBlade::alpha.size() - 1].resize(2);
    infoBlade::beta[infoBlade::beta.size() - 1].resize(2);

    std::string tempName;   
    
    readTipRadius(infoBlade::tip_radi);
    readChordLengths(infoBlade::chord);
    readInputFile(infoBlade::R, "input/degreeOfReactionConfig.dat");
    readInputFile(infoBlade::PR, "input/compressionRatioConfig.dat");
    readAngles(infoBlade::meanAlpha);
    
    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            infoBlade::diffusion[i][j] = 0.5;  

            infoBlade::alpha[i][j].resize(infoBlade::resolution);
            infoBlade::beta[i][j].resize(infoBlade::resolution);                                         
            infoBlade::lossCoefficient[i][j].resize(infoBlade::resolution);
            infoBlade::liftCoefficient[i][j].resize(infoBlade::resolution);      
            infoBlade::Mach[i][j].resize(infoBlade::resolution);
            infoBlade::pressureLoss[i][j].resize(infoBlade::resolution);         
            infoBlade::dischargeAngle[i][j].resize(infoBlade::resolution);
            infoBlade::rotateAngle[i][j].resize(infoBlade::resolution);
            infoBlade::incidenceAngle[i][j].resize(infoBlade::resolution);
            infoBlade::AoA[i][j].resize(infoBlade::resolution);
            //this doesn't work for some reasons
            //efficiency[i][j] = 0.9;
        }

        //opening a file for each batchAnalysis ofstream
        tempName = "output/batchData/" + std::to_string(i) + ".dat";
        infoBlade::batchAnalysis[i].open( tempName);
    }

    infoBlade::alpha[infoBlade::alpha.size() - 1][0].resize(infoBlade::resolution);
    infoBlade::alpha[infoBlade::alpha.size() - 1][1].resize(infoBlade::resolution);

    infoBlade::beta[infoBlade::beta.size() - 1][0].resize(infoBlade::resolution);
    infoBlade::beta[infoBlade::beta.size() - 1][1].resize(infoBlade::resolution);
    
    for(int i = 0; i < infoBlade::resolution; i++)
    {
        infoBlade::alpha[infoBlade::alpha.size() - 1][0][i] = infoBlade::meanAlpha[infoBlade::meanAlpha.size() - 1][0];
        infoBlade::alpha[infoBlade::alpha.size() - 1][1][i] = infoBlade::meanAlpha[infoBlade::meanAlpha.size() - 1][0];    
    }

    infoBlade::Area[infoBlade::Area.size()-1].resize(2);
    infoBlade::meanAlpha[infoBlade::meanAlpha.size()-1].resize(2);
    infoBlade::meanAlpha[infoBlade::meanAlpha.size()-1][0] = infoBlade::meanAlpha[infoBlade::meanAlpha.size() - 1][0];

    infoBlade::hub_radi[0][0] = infoBlade::initHub;
    infoBlade::mean_radi[0][0] = 0.5 * ( infoBlade::tip_radi[0][0] + infoBlade::hub_radi[0][0] );            

    infoBlade::Temperature[0][0] = infoBlade::T1;
    infoBlade::Pressure[0][0] = infoBlade::P1;
    infoBlade::rho[0][0] = infoBlade::rho1;

}

void getBladeElementsAngles()
{
    double dr1, dr2, radius1, radius2, y;
    double tempPhi, tempVu1, tempVu2;
    using namespace infoBlade;
    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        dr1 = ( tip_radi[i][0] - hub_radi[i][0] ) / resolution ;
        dr2 = ( tip_radi[i][1] - hub_radi[i][1] ) / resolution ;

        for(int r = 0; r <= resolution; r++)
        {
            radius1 = hub_radi[i][0] + r * dr1;
            radius2 = hub_radi[i][1] + r * dr2;
            
            y = radius1 / mean_radi[i][0]; 

            tempVu1 = (a[i] - b[i] ) /y;
            tempVu2 = (a[i] + b[i] ) /y;

            // tempVu1 = 0.0;
            // if(i < lowSize)
            // {
            // tempVu2 = work[i] / (omega1 * radius2);
            // }
            // if(i >= lowSize)
            // {
            // tempVu2 = work[i] / (omega2 * radius2);
            // }

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

            Ry[i][r] = 1.0 - ( 1.0 - R[i] ) / pow( y     , 2.0);

            //alpha[11][0][r] = meanAlpha[11][0];
            //printOut(fileOut, r, alpha[i][0][r]); 
            //printOut(fileOut2, r, alpha[i][1][r]);
            //printOut(fileOut3, r, beta[i][0][r]);
            //printOut(fileOut4, r, beta[i][1][r]);       
        }
    }
}

void thermoBlade::calculateThermoVariables()
{
    using namespace infoBlade;
    //initialising data for stage 1, the order is different from the other stages, hence the isolation
    //Initial mean radius
    mean_radi[0][0] = (tip_radi[0][0] - initHub) / 2.0 + initHub;
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

    getBladeElementsAngles();

}

void thermoBlade::calculateThermoVariables_v2()
{
    using namespace infoBlade;
    //initialising data for stage 1, the order is different from the other stages, hence the isolation
    //Initial mean radius
    mean_radi[0][0] = (tip_radi[0][0] - initHub) / 2.0 + initHub;
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

    a[0] = 0.5 * (mean_radi[0][1] + mean_radi[0][0]) * omega1 * ( 1.0 - R[0] );             
    b[0] = work[0] / ( (mean_radi[0][1] + mean_radi[0][0]) * omega1);

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

        numBlades[i][0] = 2 * PI * mean_radi[i][0] / ( chord[i][0] / solidity[i][0] );
        numBlades[i][1] = 2 * PI * mean_radi[i][1] / ( chord[i][1] / solidity[i][1] );

        chord[i][0] = solidity[i][0] * 2 * PI * mean_radi[i][0] / numBlades[i][0];
        chord[i][1] = solidity[i][1] * 2 * PI * mean_radi[i][1] / numBlades[i][1];

        //std::cout << TemperatureStag[i][0] << " " << PressureStag[i][0]  << " " << PressureStag[i][2] << std::endl;

        a[i] = 0.5 * (mean_radi[i][1] + mean_radi[i][0]) * omega1 * ( 1.0 - R[i] );             
        b[i] = work[i] / ( (mean_radi[i][1] + mean_radi[i][0]) * omega1);

        }
        //for some reasons chord[0][0] changes here
    chord[0][0] = chord[1][0];

    getBladeElementsAngles();

}

//use this  
void thermoBlade::init()
{
    sqlite3* db;
    int rc;

    resizeAndReserve();

    using namespace infoBlade;
    {
    rc = sqlite3_open("output/database/db.db", &db);
    if(rc != SQLITE_OK)
    {
        std::cout << "ERROR : failed opening the database at init;\n";
    }
    
    thermoBlade::calculateThermoVariables_v2();
    
    storeInDatabaseRecursive(db);

    sqlite3_close(db);  
    } 
}   

//only for testing
// int main()
// {       
//     infoBlade::dataBaseSetUp();
//     infoBlade::initConditionSetUp();    
//     thermoBlade::init();           

//     return 0;
// }