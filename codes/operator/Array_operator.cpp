#include "Array_operator.h"

void Array_operator::multiply(double in[2], double scalar, double out[2])
    	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = in[i]*scalar;
		}
    	};

void  Array_operator::addition(double a[2], double b[2], double out[2])
	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = a[i] + b[i];
		}
	};
	
void  Array_operator::divide(double in[2], double scalar, double out[2])
	{
		for(int i = 0; i < 2; i++)
		{
			out[i] = in[i]/scalar;
		}

	};

void  Array_operator::square(double in[2])
	{
		for(int i = 0; i < 2; i++)
		{
			in[i] = in[i] * in[i];
		}
	};
	
void  Array_operator::scalar(double in[2], double& out)
	{
		double tempX[2];
		double tempY;

		for(int i = 0; i <= 1; i++)
		{
			tempX[i] = pow(in[i], 2);
			//std::cout << i << std::endl;
			//std::cout << in[i] << " -> " << tempX[i] << std::endl;
		}
		
		tempY = tempX[0] + tempX[1];
		
		//std::cout << tempX[0] << " | " << tempX[1] << std::endl;

		out = sqrt( tempY );

		//std::cout << out << std::endl;
	};

void  Array_operator::parser(double in[2])
	{
		in[0] = sqrt(in[0]);
		in[1] = in[1]/in[0];
	};
