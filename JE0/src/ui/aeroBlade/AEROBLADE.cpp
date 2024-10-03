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
    }
}