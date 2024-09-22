#include "aeroBlade.h"

//defining all variables

std::vector<std::vector<double>> aeroBlade::I, aeroBlade::J, aeroBlade::K, aeroBlade::L;
std::vector<double> aeroBlade::Vn, aeroBlade::Vt, aeroBlade::lambda, aeroBlade::phi, aeroBlade::Sj;
std::vector<double> aeroBlade::pointX;
std::vector<double> aeroBlade::pointY;
double A, B, Cn, Dn, E;


std::complex<double> joukowskyTransform(std::complex<double> z, double shape, std::complex<double> thetaC)
{
    return exp(thetaC) * ( z + pow(shape,2) / z) ;  
}

void genBlade()
{
    std::complex<double> dummyF;
    std::complex<double> displac, aerofoil;

    double r, shape, extraR;
    r = 1.0;
    shape = 1.0;
    double multiplier = 1.0;
    double psiA, psiC, disX, disY, backFat;
    double U = 1.0;
    psiA = 0.02;
    psiC = 0.012;
    backFat = 1.02;
    disX = psiA - psiC * cos( 0 );
    disY = 0.2;
    std::ofstream output("output/misc/shape.dat");

    displac = std::complex<double>(disX, disY);
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );
    using namespace aeroBlade;
    {
        if(disX < 0)
        {
            int i = 0;
            double smallestPoint;
            for(double gh = 1.0 * ( PI ); gh >= -1.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, 0.0) / RadToDegree);
                
                pointX.insert(pointX.end(), aerofoil.real());
                pointY.insert(pointY.end(), aerofoil.imag());    

                gh -= 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }

            smallestPoint = std::min_element(pointX.begin(), pointX.end())[0];
            auto index = std::find(pointX.begin(),pointX.end(), smallestPoint);
            size_t dist = std::distance(pointX.begin(), index);

            //std::cout << dist << " " << smallestPoint << std::endl;
            //for viewing the aerofoil shape, stored in output/misc/shape.dat
            for(int i = (int)dist; i < pointX.size(); i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
            for(int i = 0; i < (int)dist; i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }

            pointX.clear();
            pointY.clear();

            std::ifstream input("output/misc/shape.dat");
            std::string temp;
            double tempD;
            for(int i = 0; i < infoBlade::resolution; i++)
            {
                std::getline(input, temp);
                std::istringstream lineStream(temp);

                lineStream >> tempD;
                pointX.insert(pointX.end(), tempD);

                lineStream >> tempD;
                pointY.insert(pointY.end(), tempD);
            }
        }
        if(disX >= 0)
        {
            int i = 0;
            double largestPoint;
            for(double gh = -2.0 * ( PI ); gh <= 0.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, 0.0) / RadToDegree);
                
                pointX.insert(pointX.end(), aerofoil.real());
                pointY.insert(pointY.end(), aerofoil.imag());           

                gh += 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }

            largestPoint = std::max_element(pointX.begin(), pointX.end())[0];
            auto index = std::find(pointX.begin(),pointX.end(), largestPoint);
            size_t dist = std::distance(pointX.begin(), index);

            //std::cout << largestPoint << std::endl;
            //for viewing the aerofoil shape, stored in output/misc/shape.dat
            for(int i = dist; i < pointX.size(); i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
            for(int i = 0; i < dist; i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }

            pointX.clear();
            pointY.clear();

            std::ifstream input("output/misc/shape.dat");
            std::string temp;
            double tempD;
            for(int i = 0; i < infoBlade::resolution; i++)
            {
                std::getline(input, temp);
                std::istringstream lineStream(temp);

                lineStream >> tempD;
                pointX.insert(pointX.end(), tempD);

                lineStream >> tempD;
                pointY.insert(pointY.end(), tempD);
            }
        }
    }
    
}

void calculateI()
{
    double temp;
    using namespace aeroBlade;
    for(int i = 0; i < infoBlade::resolution; i++)
    {
        for(int j = 0; j < infoBlade::resolution; j++)
        {
            A = -(pointX[i] - pointX[j])*cos(phi[j]) - (pointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (pointX[i] - pointX[j]), 2) + pow( (pointY[i] - pointY[j]), 2);
            Cn = sin(phi[i] - phi[j]);
            Dn = (pointY[i] - pointY[j])*cos(phi[i]) - (pointX[i] - pointX[j])*sin(phi[i]);
            E = sqrt(B - pow(A,2));
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;
            I[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
        }
    }
}

void calculateJ()
{
    double temp;
    using namespace aeroBlade;
    for(int i = 0; i < infoBlade::resolution; i++)
    {
        for(int j = 0; j < infoBlade::resolution; j++)
        {
            A = -(pointX[i] - pointX[j])*cos(phi[j]) - (pointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (pointX[i] - pointX[j]), 2) + pow( (pointY[i] - pointY[j]), 2);
            Cn = cos(phi[i] - phi[j]);
            Dn = (pointX[i] - pointX[i])*cos(phi[j]) + (pointY[i] - pointY[j])*sin(phi[i]);
            E = sqrt(B - pow(A,2));
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;
            J[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
        }
    }
}

void calculateK()
{
    double temp;
    using namespace aeroBlade;
    for(int i = 0; i < infoBlade::resolution; i++)
    {
        for(int j = 0; j < infoBlade::resolution; j++)
        {
            A = -(pointX[i] - pointX[j])*cos(phi[j]) - (pointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (pointX[i] - pointX[j]), 2) + pow( (pointY[i] - pointY[j]), 2);
            Cn = -cos(phi[i] - phi[j]);
            Dn = (pointX[i] - pointX[j])*cos(phi[i]) + (pointY[i] - pointY[j])*sin(phi[i]);
            E = sqrt(B - pow(A,2));
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;
            K[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
        }
    }
}

void calculateL()
{
    double temp;
    using namespace aeroBlade;
    for(int i = 0; i < infoBlade::resolution; i++)
    {
        for(int j = 0; j < infoBlade::resolution; j++)
        {
            A = -(pointX[i] - pointX[j])*cos(phi[j]) - (pointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (pointX[i] - pointX[j]), 2) + pow( (pointY[i] - pointY[j]), 2);
            Cn = sin(phi[i] - phi[j]);
            Dn = (pointY[i] - pointY[j])*cos(phi[i]) - (pointX[i] - pointX[j])*sin(phi[i]);
            E = sqrt(B - pow(A,2));
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;              
            J[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
        }
    }
}

//only use if disX is positive
void sourceVortexPanelMethod()
{
    using namespace aeroBlade;
    double temp;
    int i, j;
    std::vector<std::vector<double>> avg;

    std::vector<std::vector<double>> matrixA;
    std::vector<double> matrixB;

    double AoA = 0.0;
    double Vinf = 1.0;

    for(int i = 0; i < infoBlade::resolution; i++)
    {
        j = i + 1;
        if(j >= infoBlade::resolution)
        {
            j = 0;
        }
        temp = atan2( (pointY[j] - pointY[i]), (pointX[j] - pointX[i]) );
        //avg.insert(avg.end(), { pointX[i] + 0.5 * (pointX[j] - pointX[i]), pointY[i] + 0.5 * (pointY[j] - pointY[i]) });
        phi.insert(phi.end(), temp + PI / 2);       
        //std::cout << temp * RadToDegree << std::endl;
        Sj.insert(Sj.end(), sqrt( pow( pointX[j] - pointX[i], 2) + pow( pointY[j] - pointY[i], 2) ) );     
    }

    temp = atan2( (pointY[0] - pointY[infoBlade::resolution]), (pointX[0] - pointX[infoBlade::resolution]) );
    
    //the following is for checking the normal, its a sanity check
    // std::ofstream output("output/misc/shape.dat", std::ios_base::app);      
    // for(int i = 0; i < infoBlade::resolution - 1; i++)
    // {
    //     output << "\n";
    //     output << avg[i][0] << " " << avg[i][1] << std::endl;
    //     output << avg[i][0] + 0.05 * cos(phi[i]) << " " << avg[i][1] + 0.05 * sin(phi[i]) << std::endl;
    // }
    // output.close();

    calculateI();
    calculateJ();
    calculateK();
    calculateL();

    matrixA.resize(infoBlade::resolution);
    matrixB.resize(infoBlade::resolution);

    for(int i = 0; i < infoBlade::resolution; i++)
    {
        matrixA[i].resize(infoBlade::resolution);
    }

    for(int i = 0; i < infoBlade::resolution - 1; i++)
    {
        matrixB[i] = -2 * Vinf * PI * cos(phi[i] + PI / 2.0 - AoA);
         
        for(int j = 0; j < infoBlade::resolution - 1; j++)
        {
            if(i == j)
            {
                matrixA[i][i] = PI;
                continue;
            }
            matrixA[i][j] = I[i][j];
        }
    }

    double sumL;
    for(int j = 0; j < infoBlade::resolution - 1; j++)
    {
        matrixA[infoBlade::resolution - 1][j] = J[1][j] + J[infoBlade::resolution - 1][j];
        sumL += L[1][j] + L[infoBlade::resolution - 1][j];
    }

    matrixA[infoBlade::resolution - 1][infoBlade::resolution - 1] = 2 * PI - sumL;              

}   

void aeroBlade::init()
{
    
}

//only for testing
int main()
{
    infoBlade::dataBaseSetUp();
    infoBlade::initConditionSetUp();
    //thermoBlade::init();
    genBlade();
    sourceVortexPanelMethod();

    return 0;
}