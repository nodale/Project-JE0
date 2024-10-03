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
    std::ifstream readOne("input/initCondInput.dat");

    std::getline(readOne, temp);
    totalSize = std::stoi(temp);
    double tempOut[totalSize];

    std::cout << "Compression Ratios are set to : ";

    for(int i = 0; i < totalSize; i++)
    {
        stream >> tempRatio;
        if(stream.fail())
        {
            std::cout << "\nInvalid Arguments\n";
            std::cout << "Usage : CONFIGCR <num1> <num2> ... <numStages>\n";
            break;
        }
        input << tempRatio << "\n";
        std::cout << tempRatio << " ";
    }
    std::cout << "\n";

    input.close();
}

void configureInitAlpha(std::istringstream& stream)
{
    double tempAlpha;
    int totalSize;
    std::string temp;

    std::ofstream input("input/initialAlphaConfig.dat");
    std::ifstream readOne("input/initCondInput.dat");

    std::getline(readOne, temp);
    totalSize = std::stoi(temp);
    double tempOut[totalSize];

    std::cout << "Initial Alphas are set to : ";

    for(int i = 0; i <= totalSize; i++)
    {
        stream >> tempAlpha;
        if(stream.fail())
        {
            std::cout << "\nInvalid Arguments\n";
            std::cout << "Usage : CONFIGAL <num1> <num2> ... <numStages + 1>\n";
            break;
        }
        input << tempAlpha << "\n";
        std::cout << tempAlpha << " ";
    }
    std::cout << "\n";

    input.close();
}

void deleteDatabase()
{
    int rc;
    rc = std::system("rm output/database/db.db");
    
    if(rc == -1)
    {
        std::cout << "ERROR : Database has not been deleted\n";
    }
    else
    {
        std::cout << "Database has been deleted\n";
    }
}
