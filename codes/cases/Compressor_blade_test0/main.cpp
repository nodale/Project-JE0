#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <ratio>
#include <vector>

#define e 2.718281828459045
#define PI 3.14159265
#define RadToDegree 57.29577951

class Blade
{

std::ofstream fileOut;

double VX[4];
double tip_radi[9][4];
double hub_radi[9][4];
double mean_radi[9][4];
double omega1;
double omega2;
double v[9][4];
//resolution is 100 here 
double angle[9][4][100];
double Mach0;
double resolution;
double a,b;

//work-done factor
double WDF[18] = { 0.982, 0.952, 0.929, 0.910, 0.895, 0.882, 0.875, 0.868, 0.863, 0.860, 0.857, 0.855, 0.853, 0.851, 0.850, 0.849, 0.848, 0.847 };

public:

Blade(double (&tip)[9][4], double (&hub)[9][4], double w1, double res, double vx1, double vx2)
: omega1 { w1 }
, resolution { res }
, VX { vx1, vx2, vx2, vx1 }
{
    fileOut.open("out.dat");

    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            hub_radi[i][j] = hub[i][j];
            tip_radi[i][j] = tip[i][j];

            mean_radi[i][j] = ( hub[i][j] + tip[i][j] ) / 2;
        }
    }
}

/* void getExponentialVortexConst()
{
    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            auto vhub = hub_radi[i][j] * omega1;
            auto vtip = tip_radi[i][j] * omega1;

            vortexConst[i][j] = log( vtip / vhub ) / tip_radi[i][j];
            VortexVReference[i][j] = vhub;
        }
    } ctrl shift a to auto comment
} */

void getFlowPath(int i, int j, double R, double psi) //R and psi are temporarily 0.5 for the mean radius
{   
    double alpha1, alpha2, alpha3;
    double beta1, beta2, beta3;

    alpha1 = atan2( ( 1 - R - psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    alpha2 = atan2( ( 1 - R + psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    alpha3 = alpha1;

    beta1 = atan2( ( - R + psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta2 = atan2( ( - R - psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta3 = beta1;

    //std::cout << alpha1 << "\n" << alpha2 << "\n" << beta1 << "\n" << beta2 << "\n";

    getDeHallerNumber(alpha1, alpha2, beta1, beta2);
    //getSolidity( VX[j] / ( omega1 * mean_radi[i][j] ) , psi);
}

void getFlowPath(int i, int j, double R, float alpha1) //R is temporarily 0.5 for the mean radius
{   
    double phi, alpha2, alpha3;
    double beta1, beta2, beta3;

    phi = 2 * ( 1 - R - ( VX[j] * tan( alpha1 / RadToDegree ) ) / ( omega1 * mean_radi[i][j] ) );

    alpha2 = atan2( ( 1 - R + phi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    alpha3 = alpha1;

    beta1 = atan2( ( - R + phi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta2 = atan2( ( - R - phi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta3 = beta1;

    //std::cout << alpha1 << "\n" << alpha2 << "\n" << beta1 << "\n" << beta2 << "\n";

    getDeHallerNumber(alpha1, alpha2, beta1, beta2);
}

void getFlowPath(int i, int j, double R, double psi, double f) //R and psi are temporarily 0.5 for the mean radius
{   
    double alpha1, alpha2, alpha3;
    double beta1, beta2, beta3;
    double radius;

    double radi = ( tip_radi[i][j] - hub_radi[i][j] )/100;
    for(int r = 1; r <= 100; r++)
    {
    radius = ( radi * r ) + hub_radi[i][j];

    alpha1 = atan2( ( 1 - R - psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    alpha2 = atan2( ( 1 - R + psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    alpha3 = alpha1;

    beta1 = atan2( ( - R + psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta2 = atan2( ( - R - psi/2 ) , ( VX[j] / ( omega1 * mean_radi[i][j] ) ) )* RadToDegree;
    beta3 = beta1;

    //std::cout << alpha1 << "\n" << alpha2 << "\n" << beta1 << "\n" << beta2 << "\n";

    getDeHallerNumber(alpha1, alpha2, beta1, beta2);
    }
    //getSolidity( VX[j] / ( omega1 * mean_radi[i][j] ) , psi);
}
/* void getLossCoefficient(int i, int j)
{
    double lossRotor = 
} */

void getBladeAngles( int i, int j)
{
    //resolution here is temporarily 100
    double radi = ( tip_radi[i][j] - hub_radi[i][j] )/100;
    double radius, vu1, vu2, y, DegOfReaction;
   
    a = ( omega1 * mean_radi[i][j] ) * ( 0.5 - 1 );
    b = ( 0.5 * omega1 * mean_radi[i][j] * 0.5 );

    double total;

    for(int r = 1; r <= 100; r++)
    {
        radius = ( radi * r ) + hub_radi[i][j];
        y = radius / mean_radi[i][j];

        DegOfReaction = 1 + ( 0.5 * ( 1 - 2/y ) );

        if(j == 0)
        {
        vu1 = a - ( b / y );
        angle[i][j][r] = atan2( omega1 * radius - vu1, VX[j] ) * RadToDegree;
        }

        if(j == 1)
        {
        vu2 = a + ( b / y );
        angle[i][j][r] = atan2( omega1 * radius - vu2, VX[j] ) * RadToDegree;
        }

        if(j == 2)
        {
        vu2 = a + ( b / y );
        angle[i][j][r] = atan2( vu2, VX[j] ) * RadToDegree;
        }

        if(j == 3)
        {
        vu1 = a - ( b / y );
        angle[i][j][r] = atan2( vu1, VX[j] ) * RadToDegree;
        }

        total += getEnthalpyChange(radius, vu1, vu2);
        
        printOut( fileOut , radius, angle[i][j][r] );
    }
    std::cout << WDF[i]*total/100 << std::endl;

}

private:

void printOut( std::ofstream &out, double value1, double value2)
{

    out << value1 << " " << value2 << "\n";

}

void getDeHallerNumber(double alpha1, double alpha2, double beta1, double beta2)
{
    double HallerRotor = cos(alpha2 / RadToDegree )/cos(alpha1 / RadToDegree );
    double HallerStator = cos(beta1 / RadToDegree )/cos(beta2 / RadToDegree );

    std::cout << HallerRotor << " " << HallerStator << std::endl;
    /* if(HallerRotor >= 0.72 && HallerStator >= 0.72 )
    {
        std::cout << "YES" << std::endl;
    } */
    //anthing larger than 0.72 is good
}

void getSolidity(double phi, double psi)
{
    double solidity = ( 1.5 * psi ) / ( 1.55 * phi - psi );

    std::cout << solidity << "\n";
}

double getEnthalpyChange(double radius, double vu1, double vu2)
{
    return radius * omega1 * ( vu2 - vu1 );
}

double getBladeEnthalpy(double from, double to, double step, double radius, double vu1, double vu2)
{   
    //gotta fix this
    double total;
    for(int i = from; i <= to;)
    {
        total += radius * omega1 * ( vu2 - vu1 );
        i += step;
    }
    std::cout << total << std::endl;
    return total;
}

};

int main()
{
    double rTip[9][4] = {
        {
             0.08, 0.08, 0.08, 0.08
        }
    };
    double rHub[9][4] = {
        {
            0.052916782, 0.0531033701, 0.0531033701, 0.053283997
        }
    };
    double omega = 2200;
    double resol = 100;
    double v1 = 69.4377418988;
    double v2 =  v1;


    Blade test(rTip, rHub, omega, resol, v1, v2);

    test.getBladeAngles(0, 0);
    //test.getFlowPath(0, 0, 0.50, 0.50);

    return 0;
}

//plotting with gnuplot:
//plot '<datafile.dat>' with linespoints linetype 0 linewidth 2