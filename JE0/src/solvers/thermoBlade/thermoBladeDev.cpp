#include "../infoBlade/infoBlade.h"

void readTipRadius(dVec<double>& tipRadius)
{
    std::ifstream input("input/tipRadiusConfig.dat");

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            input >> tipRadius[i][j];
        }
    }
}

void readChordLengths(dVec<double>& chordLengths)
{
    std::ifstream input("input/chordConfig.dat");

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            input >> chordLengths[i][j];
        }
    }
}

void readInputFile(sVec<double>& inputMatrix, const char* argc)
{
    std::filesystem::path newpath;
    newpath = argc;
    std::ifstream input(newpath);

    std::string temp;
    int i = 0;
    while(std::getline(input, temp))
    {
        inputMatrix[i] = std::stod(temp); 
        i++;
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
    infoBlade::beta.resize(infoBlade::totalSize);
    infoBlade::meanAlpha.resize(infoBlade::totalSize + 1);
    infoBlade::meanBeta.resize(infoBlade::totalSize);
    infoBlade::PR.resize(infoBlade::totalSize);
    infoBlade::R.resize(infoBlade::totalSize);
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
        infoBlade::tip_radi[i].resize(3);
        infoBlade::tip_radi[i].resize(3);
        infoBlade::hub_radi[i].resize(3);
        infoBlade::mean_radi[i].resize(3);
        infoBlade::v[i].resize(3);
        infoBlade::diffusion[i].resize(2);
        infoBlade::numBlades[i].resize(2);

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

    infoBlade::alpha[infoBlade::alpha.size() - 1].resize(2);

    std::string tempName;   
    
    dVec<double> tipR, c;
    sVec<double> deltaP, initial_alpha1, Reaction;

    readTipRadius(tipR);
    readChordLengths(c);
    readInputFile(deltaP, "input/compressionRationConfig.dat");
    readInputFile(initial_alpha1, "input/initialAlphaConfig.dat");
    readInputFile(Reaction, "input/degreeOfReaction.dat");



    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        infoBlade::tip_radi[i] = tipR[i];
        infoBlade::PR[i] = deltaP[i];  
        infoBlade::meanAlpha[i][0] = initial_alpha1[i];    
        infoBlade::R[i] = Reaction[i];

        for(int j = 0; j < 2; j++)
        {
            infoBlade::diffusion[i][j] = 0.5;  

            infoBlade::chord[i][0] = c[i][0];      
            infoBlade::chord[i][1] = c[i][1];

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
    
    for(int i = 0; i <= infoBlade::resolution; i++)
    {
        infoBlade::alpha[infoBlade::alpha.size() - 1][0][i] = initial_alpha1[infoBlade::alpha.size() - 1];
        infoBlade::alpha[infoBlade::alpha.size() - 1][1][i] = initial_alpha1[infoBlade::alpha.size() - 1];

    }
}

void init()
{

}