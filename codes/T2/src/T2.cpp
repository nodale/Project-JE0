#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

std::ofstream outfileX("data_x.txt");
std::ofstream outfileY("data_y.txt");

double T2;
double A1 = 45.74;
double T1 = 293.15;
double Cp = 1004;
double density1 = 1.18;

int main()
{
	double n = -2.0/7.0;
	double Wi = 1000.0;
	double c1i = 1;
	double c2i = 2;
	double massflow;
	while(true)
	{
	massflow = c1i*A1*density1;
	T2 = T1 - (((Wi/massflow) - (pow(c1i,2)/2 - pow(c2i,2)/2))/Cp);
	Wi += 10;
	c1i += 0.5;
	c2i += 0.5;
	std::cout << "Power : " << Wi << "   Inlet Velocity : " << c1i << "   Outlet Velocity : " << c2i<< std::endl;
	std::cout <<  "   T2 is :" << T2 << "\n" << std::endl;

	outfileX << T2 << std::endl;
	outfileY << c1i <<std::endl;
	}	

	outfileX.close();
	outfileY.close();

	return 0;
}


