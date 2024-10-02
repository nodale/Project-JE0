#include "mainUI.h"

void quit()
{
    std::cout << "Terminating current session\n";
    std::exit(0);
}

void setupDatabase()
{
    infoBlade::dataBaseSetUp();
    std::cout << "Database is set up successfuly\n";
}

void initialiseInitCond()
{
    infoBlade::initConditionSetUp();
    std::cout << "Initial conditions are read successfuly\n";
}

void initialiseThermoVar()
{
    thermoBlade::init();
    std::cout << "Thermodynamics variables are calculated successfuly\n";
}

void initialiseThermoVar_v2()
{
    thermoBlade::init_v2();
    std::cout << "Thermodynamics variables are calculated successfuly\n";
}

void configureCompressionRatio(std::istringstream& stream)
{
    double tempRatio;
    int totalSize;
    std::string temp;

    std::ofstream input("input/compressionRatioConfig.dat");
    std::ifstream readOne("input/compressionRatioConfig.dat");

    std::getline(readOne, temp);
    totalSize = std::stoi(temp);

    for(int i = 0; i < totalSize; i++)
    {
        stream >> tempRatio;
        if(stream.fail())
        {
            std::cout << "Invalid Arguments\n";
            std::cout << "Usage : CRCONFIG <num1> <num2> ... <numN>\n";
        }
        input << tempRatio << "\n";
    }

    input.close();
}

