#include "AEROBLADE.h"
#include <fstream>

void configureBlade(std::istringstream& stream)
{
    //std::istringstream iss(args);
    double disX, disY, backFat;

    std::ofstream input("input/aerofoilConfig.dat");

    //std::cout << args << " "; 

    if(stream >> disX >> disY >> backFat)
    {
        std::cout << "The following configuration has been received : " << disX << " " << disY << " " << backFat << "\n";

        input << disX << "\n";
        input << disY << "\n";
        input << backFat << "\n";
    }
    else
    {
        std::cout << "Invalid Arguments\n";
        std::cout << "Usage : CONFIGBL <disX> <disY> <backFat>\n";
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
    }
    else
    {
        std::cout << "USAGE : CONFIRMCONFIG <numStage> <1 for rotor, 2 for stator>\n";
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

//TODO
//Maybe move getAoA somewhere else

void AEROBLADE::init()
{
    FILE* pipe = popen("gnuplot -persistent", "w");
    std::string input;
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
        if(command == "DRAWBL" or active != 0)
        {
            drawBlade(active, pipe);
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
    }
}