#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

//unit standard : m s kg

class Alpha
{
	double Area, dP;
	double  ClA, ClB[2];
	//std::vector<double> c1, c2, temp2, temp3;
	double c1[2], c2[2],temp0[2];
	double SurfaceA, SurfaceB;
	double temp1s, temp2s, temp3s, temp4s, temp5s, temp6s, temp7s;

public:

	double ClB_scalar;

	Alpha(double A, double P ,double Cl, double c_in[2], double c_out[2], double S_refA, double S_refB)
		: Area {A}
		, dP {P}
		, ClA {Cl}
		, c2 {c_in[0], c_in[1]}
		, c1 {c_out[0], c_out[1]}
		, SurfaceA {S_refA}
		, SurfaceB {S_refB}
	{

		multiply(c1, c1[0], temp0);
		scalar(temp0, temp1s);
		temp2s = Area * temp1s * (1 - dP);

		scalar(c1, temp3s);
		temp4s = 0.5 * ClA * temp3s * SurfaceA;

		scalar(c2, temp5s);
		temp6s = (-0.5) * temp5s * SurfaceB;
		
		//std::cout << temp2s << " | " << temp3s << " | " << temp5s << std::endl;

		ClB_scalar = (temp2s + temp3s)/temp6s;
		
	}

private:
	
	void multiply(double in[2], double scalar, double out[2])
	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = in[i]*scalar;
		}
	}

	void addition(double a[2], double b[2], double out[2])
	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = a[i] + b[i];
		}
	}
	
	void divide(double in[2], double scalar, double out[2])
	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = in[i]/scalar;
		}

	}

	void square(double in[2])
	{
		for(int i = 0; i < 2; i++)
		{
			in[i] = in[i] * in[i];
		}
	}
	
	void scalar(double in[2], double out)
	{
		double tempX[2];
		double tempY;

		for(int i = 0; i < 2; i++)
		{
			tempX[i] = pow(in[i], 2);
			//std::cout << in[i] << " -> " << tempX[i] << std::endl;
		}
		
		tempY = tempX[0] + tempX[1];
		
		//std::cout << tempX[0] << " | " << tempX[1] << std::endl;

		out = sqrt( tempY );
	}

	void parser(double in[2])
	{
		in[0] = sqrt(in[0]);
		in[1] = in[1]/in[0];
	}
};

int main()
{	

	std::ofstream out;	
	out.open("data/plot.dat");

	for(double var = 5; var <= 100; )
	{
	double c1_[2] = {var, var};
	double c2-[2] = {var, var};
	//the following is scaled by 10
	Alpha test_0( 0.01326,  1.05 , 0.8,  c1_,  0.0212,  0.0212 );

	out << var << " " << test_0.ClB_scalar << "\n";

	var += 0.1;
	}
	return 0;
}

//plotting with gnuplot:
//plot '<datafile.dat>' with linespoints linetype 0 linewidth 2
//
//surfaceArea = chord length for 2D
