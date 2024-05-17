#include "solver1.h"

solver1::solver1(char* file)
	: filename { file }
{
	file_out.open(file);
}

solver1::solver1(char* file, char* file2)
	: filename { file }
	, filename_real { file2 }
{
	file_out.open(file);
	file_real_out.open(file2);
}

void solver1::print_data(double x, double y, std::ofstream& target)
{
	target << x << " " << y << "\n";
}

void solver1::print_empty(std::ofstream& target)
{
	target << "\n";
}

void solver1::real(double dt, double start, double end)
{
	double tmp_y;
	for(double timer = start; timer <= end;)
	{	
  	tmp_y = func_real(timer, 0);

	print_data(timer, tmp_y, file_real_out);

	timer += dt;
	}
}	

void solver1::eulerm(double x0, double y0, double dt, double start, double end)
{
	double yb, yf, yi;
	//applying the initial conditions	
	yi = y0 + dt * func(x0, y0);
	
	//backward
	if(x0 > start)
	{
		yb = yi;

		for(double timer = (x0 - dt); timer >= start;)
		{
		
		yb = yb - dt * func(timer, yb);
		
		print_data(timer, yb, file_out);

		timer -= dt;
		}

		print_empty(file_out);
	}
	
	//forward
	
	yf = yi;
	for(double timer = (x0 + dt); timer <= end;)
	{	
		yf = yf + dt * func(timer, yf);
		
		print_data(timer, yf, file_out);

		timer += dt;
	}

	print_empty(file_out);
}

void solver1::heunm(double x0, double y0, double dt, double start, double end)
{
	double y_approx, y_approxb, y_approxf, yb, yf, yi;
	//applying the initial conditions	
	y_approx = y0 + dt * func(x0, y0);
	yi = y0 + dt * 0.5 * ( func(x0, y0) + func(x0, y_approx) );

	//backward
	if(x0 > start)
	{
		yb = yi;

		for(double timer = (x0 - dt); timer >= start;)
		{
		
		y_approxb = yb - dt * func(timer, yb);
		yb = yb - dt * 0.5 * ( func(timer, yb) + func(timer + dt, y_approxb) );

		print_data(timer, yb, file_out);

		timer -= dt;
		}

		print_empty(file_out);
	}

	//forward
	yf = yi;
	for(double timer = (x0 + dt); timer <= end;)
	{	
		y_approxf = yf + dt * func(timer, yf);
		yf = yf + dt * 0.5 * ( func(timer, yf) + func(timer + dt, y_approxf) );

		print_data(timer, yf, file_out);

		timer += dt;
	}

	print_empty(file_out);
}

uint8_t solver1::check_consistencyOrder(std::vector<double> c, std::vector<double> b, std::vector<std::vector<double>> a)
{
	uint8_t order1 = 0;
	float order2 = 0;
	double order3a = 0;
	double order3b = 0;
	for(int i = 0; i <= b.size(); i++)
	{
		order1 += b[i];
		order2 += b[i] * c[i];
		order3a += b[i] * c[i] * c[i];

		for(int j = 0; j <= b.size(); j++)
		{
			order3b += b[i] * a[i][j] * c[j];
		}
	}

	if( order1 != 1 )
	{
		return 0;
	}

	if( order2 != 0.5 )
	{
		return 1;
	}	

	if( order3a != 1/3 && order3b != 1/6 )
	{
		return 2;
	}

	return 3;

}

double solver1::ak_sum(double (&a)[4][4], double (&k)[4], int j)
{
	double sum = 0;
	for(int i = 0; i < j; i++)
	{
		sum += a[j][i] * k[i];
	}
	return sum;
}

void solver1::rungeKuttam(double x0, double y0, double dt, double start, double end, double* c, double* b, double (&a)[4][4])
{
	double sum = 0;
	double k[4] = { };
	//applying the initial conditions
	k[0] = func(x0, y0);
	k[1] = func(x0 + c[1] * dt, y0 + ( a[1][0] * k[0] ) * dt);
	k[2] = func(x0 + c[2] * dt, y0 + ( a[2][0] * k[1] + a[2][1] * k[2] ) * dt);
	k[3] = func(x0 + c[3] * dt, y0 + ( a[3][0] * k[1] + a[3][1] * k[2] + a[3][2] * k[2] ) * dt);
	double y = y0 + dt * (b[0]*k[0] + b[1]*k[1] + b[2]*k[2] + b[3]*k[3]);
	std::fill(k, k + 4, 0);

	for(double timer = x0 + dt; timer <= end;)
	
	{
		for(int j = 0; j < 4; j++)
		{
			k[j] = func(timer + c[j] * dt, y + ak_sum( a, k, j) * dt);
			sum += b[j] * k[j];
		}
/*
		k[0] = func(timer, y);
		k[1] = func(timer + c[1] * dt, y + ( a[1][0] * k[0] ) * dt);
		k[2] = func(timer + c[2] * dt, y + ( a[2][0] * k[0] + a[2][1] * k[1] ) * dt);
		k[3] = func(timer + c[3] * dt, y + ( a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2] ) * dt);
		
		sum = b[0]*k[0] + b[1]*k[1] + b[2]*k[2] + b[3]*k[3];
*/
	std::fill(k, k + 4, 0);

	y = y + dt * sum;
	print_data(timer, y, file_out);

	timer += dt;
	}

	print_empty(file_out);
}

double solver1::RK4[3][4][4] = { 
		{
			{ 0.0, 0.5 , 0.5, 1.0 },
		},

		{ 
			{ 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/3.0 },
		},

		{
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 0.5, 0.0, 0.0, 0.0 },
			{ 0.0, 0.5, 0.0, 0.0 },
			{ 0.0, 0.0, 1.0, 0.0 },
		},
};

double solver1::RK38rule[3][4][4] = { 
		{
			{ 0.0, 1.0/3.0 , 2.0/3.0, 1.0 },
		},

		{ 
			{ 1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0 },
		},

		{
			{ 0.0, 0.0, 0.0, 0.0 },
			{ 1.0/3.0, 0.0, 0.0, 0.0 },
			{ -1.0/3.0, 1.0, 0.0, 0.0 },
			{ 1.0, -1.0, 1.0, 0.0 },
		},
};
/*
std::vector<double> solver1::getRK4c = { 0, 0.5 , 0.5, 1 };
std::vector<double> solver1::getRK4b = { 1/6, 1/3, 1/3, 1/6 };
std::vector<std::vector<double>> solver1::getRK4a =
		{
			{ 0, 0, 0, 0 },
			{ 0.5, 0, 0, 0 },
			{ 0, 0.5, 0, 0 },
			{ 0, 0, 1, 0 }

		};
*/
/*template<typename d>
std::vector<std::vector<d>> solver1::linearise(std::vector<d> dFdt)
{

}*/
