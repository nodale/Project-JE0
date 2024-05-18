#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

class Array_operator
{
public :
void multiply(double in[2], double scalar, double out[2]);
void addition(double a[2], double b[2], double out[2]);
void divide(double in[2], double scalar, double out[2]);
void square(double in[2]);
void scalar(double in[2], double& out);
void parser(double in[2]);
};

