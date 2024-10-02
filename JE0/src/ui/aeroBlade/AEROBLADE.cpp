#include "AEROBLADE.h"

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
        std::cout << "Usage : CONFIG <disX> <disY> <backFat>\n";
    }

    input.close();
}

void AEROBLADE::init()
{
    std::string input;
    while(true)
    {
        std::string command, arg;
        std::vector<std::string> args;

        std::cout << "AEROBLADE> ";
        std::getline(std::cin, input);

        std::istringstream stream(input);

        stream >> command;

        if(command == "QUIT")
        {
            std::cout << "Terminating AEROBLADE\n";
            break;                                      
        }
        if(command == "CONFIG")
        {
            configureBlade(stream);
        }                               
    }
}