#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

class C2
{
	double Area;
	double ClA;
	//std::vector<double> c1, c2, temp2, temp3;
	double c1[2], temp2[2], temp3[2];
	double SurfaceA;
	double temp1;	

public:
	double c2[2];

	C2(double A, double Cl, double c_in[2], double S_ref)
		: Area {A}
		, ClA {Cl}
		, c1 {c_in[0], c_in[1]}
		, SurfaceA {S_ref}
	{
		//std::cout << Area << " " << ClA << " " << c1[0] << " " << c1[1] << " " << SurfaceA  << std::endl;

		temp1 = (ClA * SurfaceA)/(3*Area);
		
		multiply(c1, temp1, temp2);
		multiply(c1, c1[0], temp3);

		addition(temp2, temp3, c2);
		
		//std::cout << temp2[0] << " " << temp2[1] << std::endl;

		parser(c2);
	
		//std::cout << c2[0] << " " << c2[1] << std::endl;
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

	void parser(double in[2])
	{
		in[0] = sqrt(in[0]);
		in[1] = in[1]/in[0];
	}
};

int main()
{	
	//std::ofstream outx;
	//std::ofstream outy;	
	//outx.open("data/x.dat");
	//outy.open("data/y.dat");
	
	std::ofstream out;	
	out.open("data/plot.dat");

	for(double var = 0.001; var <= 0.01; )
	{
	double inlet_velocity[2] = {10, 1};

	C2 test_0(var, 0.8, inlet_velocity, 0.18);

	//outy << test_0.c2[0] <<", ";
	//outx << c1x << ", ";

	out << var << " " << test_0.c2[0] << "\n";

	var += 0.00001;
	}
	return 0;
}

//plotting with gnuplot:
//plot '<datafile.dat>' with linespoints linetype 0 linewidth 2
