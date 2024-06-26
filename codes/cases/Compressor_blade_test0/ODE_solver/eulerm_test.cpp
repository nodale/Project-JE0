#include "solver1.h"

//let x be t, and y be x
char file[] = "plot.dat";
char file_real[] = "plotreal.dat";
solver1 test1(file, file_real);

double solver1::func(double x, double y)
	{
	//	return ( (5 * x) - 3);
	//	return ( 7 * pow( y, 2) * pow( x, 3) );
		return ( 4 * pow( x, 3) * pow( 2.71828, -y) );
	};

double solver1::func_real(double x, double y)
	{ 
	//	return ( (2 * pow( 2.71828, 5 * ( x - 2 ) ))/5.0  + 3.0/5.0);
	//	return ( -1 / ( 7 * 0.25 * pow( x, 4) - 85.0/3.0 ) );
		return ( log( pow( x, 4) + pow( 2.71828, 3) - 1) );
	};

int main()
{	

	double X0 = 1;
	double Y0 = 3;

	double dt = 0.001;
	double start = -3;
	double end = 3;
	

	int co = test1.check_consistencyOrder(solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);
	std::cout << " consistency order : "<< co << std::endl;
	test1.rungeKuttam(X0, Y0, dt, start, end, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);
	
	//test1.heunrm(X0, Y0, dt, start, end);
	test1.real(0.001, start, end);

	return 0;
}

// plot '.dat' with linespoints linetype 0 linewidth 2
// plot 'plotreal.dat' with linespoints linetype 0 linewidth 2, 'plot.dat' with linespoints linetype 0 linewidth 1

