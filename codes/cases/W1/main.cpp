//TODO
//write a new runge kutta method for a vector numerical approximation

#include "solver/solver1.h"
#include "Cl_Cd/Cl_Cd.h"

#define pi 3.14159265
#define radToDeg 57.29577951
#define degToRad 0.01745329252

double c, A1, u, alpha;
double x1 = 0.0;
double x2 = 0.0;

CL_CD v( alpha, 0 );

double getBeta( double x, double y)
{
	return ( atan2(y,x) * radToDeg );
}

class X : public solver1
{
using solver1::solver1;

public:
	
	void solver1::rungeKuttam(double x0, double (&y0)[2], double dt, double start, double end, double* c, double* b, double (&a)[4][4])
{
	double y[2];
	double yi[2];
	double sum = 0;
	double k[4] = { 0 };
	//applying the initial conditions
	for(int i = 0; i < 2; i++)
	{
	k[0] = func(x0, y0[i]);
	k[1] = func(x0 + c[1] * dt, y0[i] + ( a[1][0] * k[0] ) * dt);
	k[2] = func(x0 + c[2] * dt, y0[i] + ( a[2][0] * k[1] + a[2][1] * k[2] ) * dt);
	k[3] = func(x0 + c[3] * dt, y0[i] + ( a[3][0] * k[1] + a[3][1] * k[2] + a[3][2] * k[2] ) * dt);
	yi = y0[i] + dt * (b[0]*k[0] + b[1]*k[1] + b[2]*k[2] + b[3]*k[3]);
	std::fill(k, k + 4, 0);
	
	}
	//backward
	if( x0 > start)
	{
	y[0] = yi[0];
	y[1] = yi[1];

	for(double timer = x0 - dt ; timer >= start;)
	{
		for(int i = 0; i < 2; i++)
		{
		for(int j = 0; j < 4; j++)
		{
			k[j] = func(timer + c[j] * dt, y[i] + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt);
		}
		}

	sum = b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];

	std::fill(k, k + 4, 0);

	y[0] = y[0] - dt * sum[0];
	y[1] = y[1] - dt * sum[1];
	print_data(timer, y, file_out);

	timer -= dt;
	}

	print_empty(file_out);
	}

	//forward
	y = yi;
	for(double timer = x0 +dt ; timer <= end;)
	{
		for(int j = 0; j < 4; j++)
		{
			k[j] = func(timer + c[j] * dt, y + ( a[j][0] * k[0] + a[j][1] * k[1] + a[j][2] * k[2] ) * dt);
		}

	sum = b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];

	std::fill(k, k + 4, 0);

	y = y + dt * sum;
	print_data(timer, y, file_out);

	timer += dt;
	}

	print_empty(file_out);
}
	
	//forward
	
	y = yi;
	for(double timer = (x0 + dt); timer <= end;)
	{	
		y = y + dt * func(timer, y);
		
		print_data(timer, y, file_out);

		timer += dt;
	}

	print_empty(file_out);
}
private:
	
	double func(double x, double y) final
	{
		x1 = ( 0.5 * c * x * ( v.getCl( alpha, getBeta( x , x2 ) ) * cos( getBeta( x , x2 ) * degToRad ) - v.getCd( alpha, getBeta( x , x2 ) ) * sin( getBeta( x , x2 ) * degToRad) ) );
		return x1;
	}

	double func_real(double x, double y) final
	{
		return 0;
	}
};

class Y : public solver1
{

	private:
	using solver1::solver1;

	double func(double x, double y) final
	{
		x2 = ( 0.5 * c * ( x ) * ( v.getCl( alpha, getBeta( x1 , x ) ) * sin( getBeta( x1 , x ) * degToRad) - v.getCd( alpha, getBeta( x1 , x ) ) * cos( getBeta( x1 , x ) * degToRad) ) );
		return x2;
	}

	double func_real(double x, double y) final
	{
		return 0;
	}
};


int main()
{
	char filename[] = "plot.dat";
	char filename2[] = "plot2.dat";

	X thread1(filename);
	Y thread2(filename2);

	c = 0.021213;
       	A1 = 0.013265;
       	u = 338.244809; //35.81415625;
       	alpha = 45 ;

	double X0, Y0, dt, start, end;

	X0 = 0;
	Y0 = 0;
	dt = 0.5;
	start = 0;
	end = 5;

	thread2.solver1::rungeKuttam(X0, u, dt, start, end, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);
	thread1.solver1::rungeKuttam(X0, Y0, dt, start, end, solver1::RK4[0][0], solver1::RK4[1][0], solver1::RK4[2]);
	

	return 0;
}
