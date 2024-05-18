#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdint>

class solver1
{
	double output;
	std::ofstream file_out, file_real_out;
	char* filename;
	char* filename_real;


public:
	
	solver1(char* file);
	solver1(char* file, char* file2);
	
	void real(double dt, double start, double end);

	void eulerm(double x0, double y0, double dt, double start, double end);
	void heunm(double x0, double y0, double dt, double start, double end);
	virtual void rungeKuttam(double x0, double y0, double dt, double start, double end, double* c, double* b, double (&a)[4][4]);
	
	int check_consistencyOrder(double (&c)[4], double (&b)[4], double (&a)[4][4]);

	static double RK4[3][4][4];
	static double RK38rule[3][4][4];
		
private:
	virtual double func(double x, double y);
	virtual double func_real(double x, double y);
	
	void print_data(double x, double y, std::ofstream& target);
	void print_empty(std::ofstream& target);
};
