#include <iostream>
#include <cmath>

//Cl and Cd finder for NACA65210-il

class CL_CD
{
	double Cl, Cd;
	//in degree
	double alpha, beta, gamma;
	//alpha is the angle of attack with respect to a reference coordinate
	//beta is the incoming velocity angle with respect to the reference coordinate
	//gamma is the difference between both, resulting in the angle of attack for the lift and drag calculation

public:

	CL_CD(double a, double b);
	
	double getCl(double a,double b);

	double getCd(double a,double b);

};
