#include <iomanip>
#include "solver/solver1.h"
#include "Cl_Cd/Cl_Cd.h"

#define pi 3.14159265
#define radToDeg 57.29577951
#define degToRad 0.01745329252

double c, A1, u, alpha, dt;
//k is c/2*A1 
double k;
CL_CD v( 0, 0 );

double getBeta( double x, double y)
{
	return ( atan2(y,x) * radToDeg );
}


class X : public solver1
{
double wx[2] = { 1.0 , 10.0};
double q;

public:
	
	void rungeKuttam(double x0, double y0, double dt, double start, double end, double* c, double* b, double (&a)[4][4]) final
{
	double y, yi, z, zi;
	double sum = 0;
	double sumh = 0;
	double k[4] = { 0 };
	double h[4] = { 0 };
	//applying the initial conditions

	std::fill(k, k + 4, 0);
	
	//backward
	if( x0 > start)
	{
	y = 1.0;
	z = 10;
	for(double timer = x0 ; timer >= start;)
	{
		for(int j = 0; j < 4; j++)
		{
			k[j] = func(timer + c[j] * dt, y + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt, z + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2]))[0];
			h[j] = func(timer + c[j] * dt, y + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt, z + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2]))[1];
		}

	sum = b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];
	sumh = b[0] * h[0] + b[1] * h[1] + b[2] * h[2] + b[3] * h[3];
	std::fill(k, k + 4, 0);

	y = y - dt * sum;
	z = z - dt * sumh;
	print_data(timer, y, file_out);

	timer -= dt;
	}

	print_empty(file_out);
	}

	//forward
	y = 1.0;
	z = 10;
	for(double timer = x0  ; timer <= end;)
	{
		for(int j = 0; j < 4; j++)
		{
			k[j] = func(timer + c[j] * dt, y + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt, z + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2]))[0];
			h[j] = func(timer + c[j] * dt, y + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt, z + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2]))[1];
		}

	sum = b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];
	sumh = b[0] * h[0] + b[1] * h[1] + b[2] * h[2] + b[3] * h[3];
	std::fill(k, k + 4, 0);

	y = y + dt * sum;
	z = z + dt * sumh;
	std::cout << y << " " << z << " " << std::endl;
	print_data(timer, y, file_out);

	timer += dt;
	}

	print_empty(file_out);
}

private:
	using solver1::solver1;
	
	/*
	double func(double x, double y) final
	{
		q = k * sqrt( pow( y, 2 ) + pow( wx[1], 2 ) ) / wx[0];

		wx[0] = q * ( v.getCl( alpha, getBeta( y , wx[1] ) ) - v.getCd( alpha, getBeta( y , wx[1] ) ) );
		return wx[0];
	}
	*/

	
	double* func(double x, double y, double z) final
	{
		//y is x1, z is x2	
		q = k * sqrt( pow( y, 2 ) + pow( z, 2 ) ) / y;

		wx[0] = q * ( v.getCl( alpha, getBeta( y , z ) ) - v.getCd( alpha, getBeta( y , z ) ) );
		wx[1] = q * ( -v.getCl( alpha, getBeta( y , z ) ) - v.getCd( alpha, getBeta( y , z ) ) );
		return wx;
	}

/*
 		 q = k * sqrt( pow( y, 2 ) + pow( x2, 2 ) ) / x1;

		x1 = q * ( v.getCl( alpha, getBeta( y , x2 ) ) - v.getCd( alpha, getBeta( y , x2 ) ) );
		x2 = q * ( v.getCl( alpha, getBeta( y , x2 ) ) + v.getCd( alpha, getBeta( y , x2 ) ) );
 */
	double func_real(double x, double y) final
	{
		return 0;
	}
};

/*
class Y : public solver1
{

	private:
	using solver1::solver1;

	double func(double x, double y) final
	{
		
		return 0;
	}

	double func_real(double x, double y) final
	{
		return 0;
	}
};
*/

int main()
{
	std::cout << std::fixed << std::setprecision(5);

	char filename[] = "plot.dat";
	//char filename2[] = "plot2.dat";

	X thread1(filename);
	//Y thread2(filename2);

	c = 0.021213;
       	A1 = 0.013265;
       	u = 338.244809; //35.81415625;
       	alpha = 45.0 ;

	double X0, Y0, start, end;

	X0 = 0.0;
	Y0 = 0.0;
	dt = 0.01;
	start = 0.0;
	end = 1000.0;
	k = c / ( 2.0 * A1);

	//thread2.solver1::rungeKuttam(X0, u, dt, start, end, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);
	thread1.solver1::rungeKuttam(X0, Y0, dt, start, end, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);

	

	return 0;
}
