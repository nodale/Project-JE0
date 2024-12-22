#include "AEROBLADE.h"
#include <fstream>

void configureBlade(std::istringstream& stream)
{
    double disX, disY, backFat;

    std::ofstream input("input/aerofoilConfig.dat");

    if(stream >> disX >> disY >> backFat)
    {
        std::cout << "The following configuration has been received : " << disX << " " << disY << " " << backFat << "\n";

        input << disX << "\n";
        input << disY << "\n";
        input << backFat << "\n";
    }
    else
    {   
        stream.clear();

        std::cout << "Invalid Arguments\n";
        std::cout << "USAGE : CONFIGBL <disX> <disY> <backFat>\n";
    }

    input.close();
}

void drawBlade(int active, FILE* pipe)
{
    aeroBlade::genBlade(0,0);

    if(active == 0)
    {
    fprintf(pipe, "set title 'Aerofoil Profile'\n");
    fprintf(pipe, "set xrange[-0.5:1.5]\n");
    fprintf(pipe, "set yrange[-0.5:0.5]\n");
    fprintf(pipe, "plot 'output/misc/shape.dat' with linespoints linetype -1 linewidth 2\n");
    fflush(pipe); 
    fprintf(pipe, "\n");
    //pclose(pipe);
    }
    if(active == 1)
    {
        fprintf(pipe, "set xrange[-0.5:1.5]\n");
        fprintf(pipe, "set yrange[-0.5:0.5]\n");
        fprintf(pipe, "replot\n");
        fflush(pipe);
    }
    if(active == 2)
    {
        fprintf(pipe, "quit\n");
        fflush(pipe);
        pclose(pipe);
    }
}

void storeConfig(std::istringstream& stream)
{
    double disX, disY, backFat;

    int stage, j;

    std::string temp1;

    std::ifstream input("input/aerofoilConfig.dat");

    std::getline(input, temp1);
    disX = std::stod(temp1);
    std::getline(input, temp1);

    if(j == 1)
    {
        disY = -std::stod(temp1);
    }
    if(j == 0)
    {
        disY = std::stod(temp1);                                
    }

    std::getline(input, temp1);
    backFat = std::stod(temp1);

    if(stream >> stage >> j)
    {
        sqlite3* db;
        sqlite3_open("output/database/db.db", &db);

        std::string text, rotorOrStator;

        if(j == 0)
        {
            rotorOrStator = "_rotor";
        }
        if(j == 1)
        {
            rotorOrStator = "_stator";
        }

        text = "disX" + rotorOrStator;
        infoBlade::storeInDesignDatabase(db, text, disX, stage);

        text = "disY" + rotorOrStator;
        infoBlade::storeInDesignDatabase(db, text, disY, stage);

        text = "backFat" + rotorOrStator;
        infoBlade::storeInDesignDatabase(db, text, backFat, stage);
    
        std::cout << "Successfully stored the blade configuration for stage number " << stage << std::endl; 
    }
    else
    {
        std::cout << "USAGE : CONFIRMCONFIG <numStage> <0 for rotor, 1 for stator>\n";
    }
}

void findRandomCombinationAlpha(std::istringstream& stream)
{
    int sampleSize, maxTries;

    if(stream >> sampleSize >> maxTries)
    {
        std::cout << "Finding the best random alpha1 combination out of " << sampleSize << " samples\n";
        aeroBlade::findCombinationAlpha(sampleSize,maxTries);
    }
    else
    {
        std::cout << "USAGE : RANDOMALPHA <sampleSize> <maxAttempts>\n";
    }
}

void findRandomCombinationFull(std::istringstream& stream)
{
    int sampleSize, maxTries;

    if(stream >> sampleSize >> maxTries)
    {
        std::cout << "Finding the best random alpha1 and omega1 combination out of " << sampleSize << " samples\n";
        aeroBlade::findCombinationFull(sampleSize,maxTries);
    }
    else
    {
        std::cout << "USAGE : RANDOMALPHA <sampleSize> <maxAttempts>\n";
    }
}

void drawLiftCoefficients(std::istringstream& stream)
{
    FILE* pipe = popen("gnuplot -persistent", "w");

    fprintf(pipe, "set title 'Lift Coefficients'\n");
    fprintf(pipe, "set xrange[-0.2:1.2]\n");
    fprintf(pipe, "set yrange[-2.0:2.0]\n");
    fprintf(pipe, "plot 'plot' with linespoints linetype -1 linewidth 2\n");

    for(int j = 0; j < 2; j++)
    {
        for(int i = 0; i < infoBlade::totalSize; i++)
        {
            for(int r = 0; r < infoBlade::resolution; r++)
            {
                fprintf(pipe, "%lf %lf\n", (double)r/infoBlade::resolution,infoBlade::liftCoefficient[i][j][r]);
            }
            fprintf(pipe, "\n");
        }
        fprintf(pipe, "\n");
    }

    fflush(pipe); 
    pclose(pipe);
}

void drawDeHallers(std::istringstream& stream)
{
    FILE* pipe = popen("gnuplot -persistent", "w");

    fprintf(pipe, "set title 'De Hallers Number'\n");
    fprintf(pipe, "set xrange[-0.1:1.1]\n");
    fprintf(pipe, "set yrange[-0.0:2.0]\n");
    fprintf(pipe, "plot '-' with linespoints linetype -1 linewidth 2\n");

    double theta = 0.0;

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int r = 0; r < infoBlade::resolution; r++)
        {
            theta = cos(infoBlade::beta[i][0][r]/RadToDegree) / cos(infoBlade::beta[i][1][r]/RadToDegree); 
            fprintf(pipe, "%lf %lf\n", (double)r/infoBlade::resolution, theta);
        }
        fprintf(pipe, "\n");
    }

    fprintf(pipe, "\n");

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int r = 0; r < infoBlade::resolution; r++)
        {
            theta = cos(infoBlade::alpha[i+1][0][r]/RadToDegree) / cos(infoBlade::alpha[i][1][r]/RadToDegree); 
            fprintf(pipe, "%lf %lf\n", (double)r/infoBlade::resolution, theta);
        }
        fprintf(pipe, "\n");
    }
    fprintf(pipe, "\n");

    fflush(pipe); 
    pclose(pipe);
}

void drawAngles(std::istringstream& stream)
{
    FILE* pipe = popen("gnuplot -persistent", "w");

    fprintf(pipe, "set title 'Blade Elements Angle'\n");
    fprintf(pipe, "set xrange[-0.1:1.1]\n");
    fprintf(pipe, "set yrange[-100.0:100.0]\n");
    fprintf(pipe, "plot '-' with linespoints linetype -1 linewidth 2\n");

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int r = 0; r < infoBlade::resolution; r++)
        {
            fprintf(pipe, "%lf %lf\n", (double)r/infoBlade::resolution, infoBlade::alpha[i][0][r]);
        }
        fprintf(pipe, "\n");
    }

    fprintf(pipe, "\n");

    for(int i = 0; i < infoBlade::totalSize; i++)
    {
        for(int r = 0; r < infoBlade::resolution; r++)
        {
            fprintf(pipe, "%lf %lf\n", (double)r/infoBlade::resolution, infoBlade::alpha[i][1][r]);
        }
        fprintf(pipe, "\n");
    }
    fprintf(pipe, "\n");

    fflush(pipe); 
    pclose(pipe);
}

void runSim(std::istringstream& stream)
{
    int i, j;

    if(stream >> i >> j)
    {
        simBlade::generateAerofoilModel( i - 1, j);

        system("cp -r output/systemFile/* simCase/system\n");
    }
    else        
    {
        std::cout << "USAGE : RUNSIM <numStage> < 0 for rotor; 1 for stator >\n";
    }
}

//TODOs
//code functions for running CFD

void AEROBLADE::init()
{
    FILE* pipe = popen("gnuplot -persistent", "w");
    std::string input;

    bool initiate = 0;
    
    while(true)
    {
        std::string command, arg;
        std::vector<std::string> args;

        std::cout << "AEROBLADE> ";
        std::getline(std::cin, input);

        std::istringstream stream(input);

        int active = 0;

        stream >> command;

        if(command == "QUIT")
        {
            std::cout << "Terminating AEROBLADE\n";
            break;                                      
        }
        if(command == "STOPDRAWING")
        {
            active = 2;
        }
        if(command == "CONFIGBL")
        {
            configureBlade(stream);
            active = 1;
        }   
        if(command == "DRAWBL" or (active != 0 && initiate == 1))
        {
            drawBlade(active, pipe);
            initiate = 1;
        }           
        if(command == "CONFIRMCONFIG")
        {
            storeConfig(stream);
        }    
        if(command == "RANDOMCOMBALPHA")
        {   
            findRandomCombinationAlpha(stream);
        }
        if(command == "RANDOMCOMBFULL")
        {   
            findRandomCombinationFull(stream);
        }  
        if(command == "DRAWCL")
        {
            drawLiftCoefficients(stream);
        }         
        if(command == "DRAWDH")
        {
            drawDeHallers(stream);
        }
        if(command == "DRAWALPHAS")
        {
            drawAngles(stream); 
        }
        if(command == "RUNSIM")
        {
            runSim(stream);
        }
    }
}