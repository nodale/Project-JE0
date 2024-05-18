#include <iostream>
#include <fstream>
#include "operator/Array_operator.h"
#include "Cl_Cd/Cl_Cd.h"

#define PI 3.1416

class Stator : public Array_operator, CL_CD
{
	double Area, dP, SurfaceA, Cl, Cd;
	double c2[2], c3[2];
	
	std::ofstream out;

public:

	template <typename d>
	Compressor( d area, d p, d clA, double c_in[2], d s_refA )
		: Area {area}
		, dP {p}
		, c2 {c_in[0], c_in[1]}
		, SurfaceA {s_refA}
	{


	}
	
private:
	
	template <typename d>
	void find_c3(d Area, d dP, double c2[2], d SurfaceA)
	{
		d tmp0 = (0.5 * SurfaceA)/(dP * Area);
		d tmp1 = sqrt( pow(c2[0], 2) + pow(c2[1], 2) );

		d tmp2 = atan2(c2[1]/c2[0]);
		CL_CD param(45, tmp2);
		
		d tmp3 = (param.Cl * cos( tmp2 * PI / 180.0)) - (param.Cd * sin( tmp2 * PI / 180.0));a

		d tmp4 = (tmp0 * tmp1 * tmp3) + (c2[0] * c2[1]);
	}
	
};

int main()
{
	std::ofstream out;
	out.open("data/plot.dat");

	double area = 0.01326;
	double p_gain = 1.05;
	double clA = 0.8;
	double s_refA = 0.0212;
	double s_refB = 0.0212;

	for(double var = 10; var <= 200;)
	{
	
	double c1[2] = {var, var};

	Stator Stage1(area, p_gain, clA, c1, s_refA, s_refB);
	
	out << var << " " << Stage1.c2[0] << " "<< Stage1.ClB << "\n";

	var += 0.1;
	}

	return 0;
}
