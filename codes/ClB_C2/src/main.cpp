#include <iostream>
#include <fstream>
#include "operator/Array_operator.h"

class Compressor : public Array_operator
{
	double Area, dP, ClA, SurfaceA, SurfaceB;
	double c1[2], temp0[2];
	
	std::ofstream out;

public:
	
	double ClB;
	double c2[2];

	template <typename d>
	Compressor( d area, d p, d clA, double c_in[2], d s_refA, d s_refB )
		: Area {area}
		, dP {p}
		, ClA {clA}
		, c1 {c_in[0], c_in[1]}
		, SurfaceA {s_refA}
		, SurfaceB {s_refB}
	{

		find_c2(Area, ClA, c1, SurfaceA);
		find_ClB(Area, dP, ClA, c1, c2, SurfaceA, SurfaceB);

	}
	
private:
	
	template <typename d>
	void find_c2(d Area , d ClA, double c1[2], d SurfaceA)
	{
		d temp1 = (ClA * SurfaceA)/(3 * Area);
		
		d temp2[2], temp3[2];
		multiply(c1, temp1, temp2);
		multiply(c1, c1[0], temp3);
		
		addition(temp2, temp3, c2);

		parser(c2);
	}
	
	template <typename d>
	void find_ClB(d Area, d dP, d ClA, double c1[2], double c2[2], d SurfaceA, d SurfaceB)
	{
		d temp0[2], temp1, temp2;
		scalar(c1, temp1);
		temp2 =  Area * temp1 * (1 - dP) * c1[0];
		
		d temp3, temp4;
		scalar(c1, temp3);
		
		temp4 = 0.5 * ClA * temp3 * SurfaceA;
						
		d temp5, temp6;
		scalar(c2, temp5);
		temp6 = (-0.5) * temp5 * SurfaceB;

		ClB = (temp2 - temp3)/temp6;

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

	Compressor Stage1(area, p_gain, clA, c1, s_refA, s_refB);
	
	out << var << " " << Stage1.c2[0] << " "<< Stage1.ClB << "\n";

	var += 0.1;
	}

	return 0;
}
