#include "main/mainUI.h"
#include "aeroBlade/AEROBLADE.h"

int main()
{
    std::string input;

    while(true)
    {
        std::string command, arg;
        std::vector<std::string> args;

        std::cout << "JE0> ";
        std::getline(std::cin, input);

        std::istringstream stream(input);

        stream >> command;

        if(command == "QUIT")
        {
            quit();
        }
        if(command == "SETUP")
        {
            setupDatabase();
        }
        if(command == "INIT")
        {
            initialiseInitCond();
        }
        if(command == "THERMO")
        {
            initialiseThermoVar();
        }
        if(command == "THERMO_V2")
        {
            initialiseThermoVar_v2();
        }
        if(command == "AEROBLADE")
        {
            AEROBLADE::init();
        }
        if(command == "CONFIGCR")
        {
            configureCompressionRatio(stream);
        }
        // else
        // {
        //     std::cout << "Unknown command\n";
        // }
        
        
    }

    return 0;
}