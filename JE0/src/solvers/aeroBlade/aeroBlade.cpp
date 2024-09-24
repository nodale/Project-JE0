#include "aeroBlade.h"
#include <cmath>

//defining all variables

dVec<double> aeroBlade::I, aeroBlade::J, aeroBlade::K, aeroBlade::L;
sVec<double> aeroBlade::Vn, aeroBlade::Vt, aeroBlade::lambda, aeroBlade::phi, aeroBlade::Sj;
sVec<double> aeroBlade::pointX;
sVec<double> aeroBlade::pointY;
sVec<double> aeroBlade::midPointX;
sVec<double> aeroBlade::midPointY;
double aeroBlade::aeroGamma;
double A, B, Cn, Dn, Ct, Dt, E;

//function to inverse matrix A by means of LU decomposition method
dVec<double> inverseLU(const dVec<double>& A)
{
    std::vector<std::vector<double>> decU, decL, inverseL, inverseU, output;
    double temp = 0;
    decU.resize(A.size());
    decL.resize(A.size());
    inverseL.resize(A.size());
    inverseU.resize(A.size());
    output.resize(A.size());

    for(int i = 0; i < decU.size(); i++)
    {
        decU[i].resize(A.size());
        decL[i].resize(A.size());
        inverseL[i].resize(A.size());
        inverseU[i].resize(A.size());
        output[i].resize(A.size());
    }

    //setting up Dooltille's algorithm                      
    for(int i = 0; i < decL.size(); i++)
    {
        decL[i][i] = 1.0;
        decU[0][i] = A[0][i];  
    }

    for(int j = 0; j < decU.size(); j++)
    {
        for(int i = 0; i <= j; i++)
        {
            for(int k = 0; k <= i - 1; k++)
            {
                temp += decL[i][k] * decU[k][j];
            }
            decU[i][j] = A[i][j] - temp;
            temp = 0;
        }
    }

    for(int i = 0; i < decL.size(); i++)
    {
        for(int j = 0; i > j; j++)
        {
            if(i == j)
            {
                continue;
            }
            for(int k = 0; k <= j - 1; k++)
            {
                temp += decL[i][k] * decU[k][j];                
            }
            decL[i][j] = ( 1.0 / decU[j][j]) * ( A[i][j] - temp );
            temp = 0;
        }
    }

    for(int j = 0; j < decU.size(); j++)
    {
        for(int k = 0; k < decU.size() - 1; k++)
        {
            temp += decL[j][k] * decU[k][j];
        }
        decU[j][j] = A[j][j] - temp;
        temp = 0;
    }

    for(int i = 0; i < inverseL.size(); i++)
    {
        inverseL[i][i] = 1.0 / decL[i][i];
        inverseL[i][0] = (-1.0 / decL[i][i]) * decL[i][0] * inverseL[0][0]; 
        inverseU[i][i] = 1.0 / decU[i][i];
    }
    for(int i = 0; i < inverseU.size(); i++)
    {
        for(int n = inverseU.size() - 1; n >= 0; n--)
        {   
            for(int k = i + 1; k <= n; k++)
            {
                temp += decU[i][k] * inverseU[k][n];    
            }
            inverseU[i][n] = (-1.0 / decU[i][i]) * temp;
            temp = 0;
        }
    }
    for(int i = 0; i < inverseL.size(); i++)
    {
        for(int j = 0; i > j; j++)
        {
            for(int k = j; k <= i - 1; k++)
            {
                temp += decL[i][k] * inverseL[k][j];
            }
            inverseL[i][j] = -temp / decL[j][j];
            temp = 0;
        }
    }
    for(int j = 0; j < inverseU.size(); j++)
    {
        for(int i = 0; i <= j; i++)
        {
            if(i == j)
            {
                continue;
            }
            for(int k = i + 1; k <= j; k++)
            {
                temp += decU[i][k] * inverseU[k][j];
            }
            inverseU[i][j] = -temp / decU[i][i];
            temp = 0;
        }
    }
    std::ofstream out("output/misc/inverseA.dat");
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A.size(); j++)
        {
            for(int k = 0; k < A.size(); k++)
            {
                temp += inverseU[i][k] * inverseL[k][j];
            }
            output[i][j] = temp;
            out << output[i][j] << " ";
            temp = 0;
        }
        out << "\n";
    }

    return output;
}

dVec<double> inverseLU_v2(const dVec<double>& A)
{
    int n = A.size();
    dVec<double> decU(n, sVec<double>(n, 0));
    dVec<double> decL(n, sVec<double>(n, 0));
    dVec<double> inverseL(n, sVec<double>(n, 0));
    dVec<double> inverseU(n, sVec<double>(n, 0));
    dVec<double> output(n, sVec<double>(n, 0));

    for (int i = 0; i < n; i++) 
    {
        decL[i][i] = 1.0; 

        for (int j = 0; j < n; j++) 
        {
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += decL[i][k] * decU[k][j];
            }
            decU[i][j] = A[i][j] - sum;    

            if (i > j) 
            {
                sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += decL[i][k] * decU[k][j];
                }
                decL[i][j] = (A[i][j] - sum) / decU[j][j];
            }
        }
    }

    //inverts L matrix using forward substitution
    for (int i = 0; i < n; i++) 
    {
        inverseL[i][i] = 1.0;
        for (int j = 0; j < i; j++) 
        {
            double sum = 0;
            for (int k = j; k < i; k++) 
            {
                sum += decL[i][k] * inverseL[k][j];
            }
            inverseL[i][j] = -sum;
        }
    }

    //inverts U matrix using backward substitution
    for (int i = n - 1; i >= 0; i--) 
    {
        inverseU[i][i] = 1.0 / decU[i][i];
        for (int j = i + 1; j < n; j++) 
        {
            double sum = 0;
            for (int k = i + 1; k <= j; k++) 
            {
                sum += decU[i][k] * inverseU[k][j];
            }
            inverseU[i][j] = -sum / decU[i][i];
        }
    }

    //outstream for double check
    std::ofstream out("output/misc/inverseA.dat");
    //mltiplies inverseU and inverseL to get inverse of A
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            double sum = 0;
            for (int k = 0; k < n; k++) 
            {
                sum += inverseU[i][k] * inverseL[k][j];
            }
            output[i][j] = sum;
            out << output[i][j] << " ";
        }
        out << "\n";
    }

    return output;
}

//function to multiply matrix inverse A with B
sVec<double> multiplyAB(const dVec<double>& invA, const sVec<double>& B)
{
    sVec<double> C;

    double temp = 0;
    for(int i = 0; i < invA.size() - 1; i++)
    {
        for(int j = 0; j < B.size(); j++)
        {
            temp += invA[i][j] * B[j];
        }
        C.insert(C.end(), temp);
        temp = 0;
    }

    int i = invA.size() - 1;    
    for(int j = 0; j < B.size(); j++)
    {
        temp += invA[i][j] * B[j];
    }
    aeroBlade::aeroGamma = temp;
    temp = 0;

    double tempSum = 0;
    for(int i = 0; i < C.size(); i++)
    {
        tempSum += C[i];
    }
    std::cout << tempSum << std::endl;
    std::cout << aeroBlade::aeroGamma << std::endl;
    return C;
}

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
    disY = 0.034;
    std::ofstream output("output/misc/shape.dat");
                
    aeroBlade::pointX.resize(infoBlade::resolution);
    aeroBlade::pointY.resize(infoBlade::resolution);

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
                
                pointX[i] = aerofoil.real();
                pointY[i] = aerofoil.imag();    

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
                
                pointX[i] = aerofoil.real();
                pointY[i] = aerofoil.imag();              

                gh += 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }

            largestPoint = std::max_element(pointX.begin(), pointX.end())[0];
            auto index = std::find(pointX.begin(),pointX.end(), largestPoint);
            size_t dist = std::distance(pointX.begin(), index);

            //std::cout << largestPoint << std::endl;
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
        }
    }
    
}

void calculateIJ()
{
    //std::ofstream out("output/misc/I.dat");
    double temp, temp2;
    using namespace aeroBlade;
    for(int i = 0; i < I.size(); i++)
    {
        for(int j = 0; j < I.size(); j++)
        {
            A = -(midPointX[i] - pointX[j])*cos(phi[j]) - (midPointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (midPointX[i] - pointX[j]), 2) + pow( (midPointY[i] - pointY[j]), 2);
            Cn = sin(phi[i] - phi[j]);
            Dn = (midPointY[i] - pointY[j])*cos(phi[i]) - (midPointX[i] - pointX[j])*sin(phi[i]);
            Ct = -cos(phi[i] - phi[j]);
            Dt = (midPointX[i] - pointX[j])*cos(phi[i]) + (midPointY[i] - pointY[j])*sin(phi[i]);
            E = sqrt( B - pow(A,2) );
            if( (B - pow(A,2)) < 0 or std::isnan(E) == 1 or std::isinf(E) == 1)
            {
                E = 0;  
            }
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;
            I[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
            J[i][j] = 0.5 * Ct * log( temp / B ) + ( Dt - A * Ct ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;

            if(std::isnan(I[i][j]) == 1 || std::isinf(I[i][j]) == 1)
            {
                I[i][j] = 0;
            }
            if(std::isnan(J[i][j]) == 1 || std::isinf(J[i][j]) == 1)
            {
                J[i][j] = 0;
            }   
            if(i == j)
            {
                I[i][j] = 0.0;
                J[i][j] = 0.0;
            }

            //out << I[i][j] << " ";
        }
        //out << "\n";
    }
}

void calculateKL()
{
    //std::ofstream out("output/misc/K.dat");

    double temp, temp2;
    using namespace aeroBlade;
    for(int i = 0; i < K.size(); i++)  
    {
        for(int j = 0; j < K.size(); j++)
        {
            A = -(midPointX[i] - pointX[j])*cos(phi[j]) - (midPointY[i] - pointY[j])*sin(phi[j]);
            B = pow( (midPointX[i] - pointX[j]), 2) + pow( (midPointY[i] - pointY[j]), 2);
            Cn = -cos(phi[i] - phi[j]);
            Dn = (midPointX[i] - pointX[j])*cos(phi[i]) + (midPointY[i] - pointY[j])*sin(phi[i]);
            Ct = sin(phi[j] - phi[i]);
            Dt = -(midPointY[i] - pointY[j])*cos(phi[i]) + (midPointX[i] - pointX[j])*sin(phi[i]);
            E = sqrt(  std::max( B - pow(A,2), 0.0 ) );
            if( (B - pow(A,2)) < 0 )
            {
                E = 0;
            }
            temp = pow( Sj[j], 2) + 2 * A * Sj[j] + B;
            K[i][j] = 0.5 * Cn * log( temp / B ) + ( Dn - A * Cn ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;
            L[i][j] = 0.5 * Ct * log( temp / B ) + ( Dt - A * Ct ) * ( atan2( ( Sj[j] + A ) , E ) - atan2( A , E ) ) / E;

            if(std::isnan(K[i][j]) == 1 || std::isinf(K[i][j]) == 1)
            {
                K[i][j] = 0;
            }
            if(std::isnan(L[i][j]) == 1 || std::isinf(L[i][j]) == 1)
            {
                L[i][j] = 0;
            }
            if(i == j)
            {
                K[i][j] = 0.0;
                L[i][j] = 0.0;
            }
            //out << K[i][j] << " ";
        }
        //out << "\n";
    }
}

//only use if disX is positive
void sourceVortexPanelMethod()
{
    using namespace aeroBlade;
    double temp;
    dVec<double> avg;
    dVec<double> matrixA, inverseMatrixA;
    sVec<double> matrixB, aeroCp, aeroCN, aeroCA, aeroCm;

    double AoA = 0.0;
    double Vinf = 1.0;

    std::ifstream input("output/misc/shape.dat");
    std::string tempS;
    double tempD;
    int numPanels = infoBlade::resolution - 1;
    int numBoundaries = infoBlade::resolution;
    for(int i = 0; i < numBoundaries; i++)                                                  
    {
        std::getline(input, tempS);
        std::istringstream lineStream(tempS);

        lineStream >> tempD;
        pointX.insert(pointX.end(), tempD);

        lineStream >> tempD;
        pointY.insert(pointY.end(), tempD);
    }

    for(int i = 0; i < numPanels; i++)
    {
        int j = i + 1;

        temp = atan2( (pointY[j] - pointY[i]), (pointX[j] - pointX[i]) );
        //avg.insert(avg.end(), { pointX[i] + 0.5 * (pointX[j] - pointX[i]), pointY[i] + 0.5 * (pointY[j] - pointY[i]) });
        midPointX.insert(midPointX.end(), pointX[i] + 0.5 * (pointX[j] - pointX[i]) );
        midPointY.insert(midPointY.end(), pointY[i] + 0.5 * (pointY[j] - pointY[i]) );
        phi.insert(phi.end(), temp);       
        //std::cout << temp * RadToDegree << std::endl;
        Sj.insert(Sj.end(), sqrt( pow( pointX[j] - pointX[i], 2) + pow( pointY[j] - pointY[i], 2) ) );     
    }
    //temp = atan2( (pointY[0] - pointY[infoBlade::resolution - 1]), (pointX[0] - pointX[infoBlade::resolution - 1]) );
    
    //the following is for checking the normal, its a sanity check
    // std::ofstream output("output/misc/shape.dat", std::ios_base::app);      
    // for(int i = 0; i < infoBlade::resolution - 1; i++)
    // {
    //     output << "\n";
    //     output << avg[i][0] << " " << avg[i][1] << std::endl;
    //     output << avg[i][0] + 0.05 * cos(phi[i]) << " " << avg[i][1] + 0.05 * sin(phi[i]) << std::endl;
    // }
    // output.close();

    I.resize(numPanels);
    J.resize(numPanels);
    K.resize(numPanels);
    L.resize(numPanels);
    lambda.resize(numPanels);
    for(int i = 0; i < I.size(); i++)
    {
        I[i].resize(numPanels);
        J[i].resize(numPanels);
        K[i].resize(numPanels);
        L[i].resize(numPanels);
    }

    calculateIJ();   
    calculateKL();

    matrixA.resize(numPanels + 1);
    matrixB.resize(numPanels + 1);

    for(int i = 0; i < numPanels + 1; i++)
    {
        matrixA[i].resize(numPanels + 1);
    }

    for(int i = 0; i < numPanels; i++)
    {
        matrixB[i] = -2 * Vinf * PI * cos(phi[i] + PI / 2.0 - AoA);
        for(int j = 0; j < numPanels; j++)
        {
            matrixA[i][j] = I[i][j];
            matrixA[i][numPanels] += -K[i][j];  
            //std::cout << matrixA[i][j] << " " << i << " " << j << std::endl;      
        }
        matrixA[i][i] = PI;
    }

    matrixB[numPanels] = -Vinf * 2 * PI * ( sin(phi[0] + PI / 2.0 - AoA) + sin(phi[numPanels - 1] + PI / 2.0 - AoA));

    double sumL = 0;
    for(int j = 0; j < numPanels; j++)
    {
        sumL += L[0][j] + L[numPanels - 1][j];
        matrixA[numPanels][j] = J[0][j] + J[numPanels - 1][j];
    }

    matrixA[numPanels][numPanels] = 2 * PI - sumL;        
    inverseMatrixA = inverseLU_v2(matrixA);     

    //outstream for sanity check
    std::ofstream out2("output/misc/matrixA.dat");
    for(int i = 0; i <= numPanels; i++)
    {
        for(int j = 0; j <= numPanels; j++)
        {
            out2 << matrixA[i][j] << " ";
        }
        out2 << "\n";
    }

    aeroBlade::lambda = multiplyAB(inverseMatrixA, matrixB);
    double sumV = 0;
    for(int i = 0; i < numPanels; i++)
    {
        for(int j = 0; j < numPanels; j++)
        {
            sumV += lambda[j] * J[i][j] / ( 2.0 * PI ) - aeroBlade::aeroGamma * L[i][j] / ( 2.0 * PI );
        }
        Vt.insert(Vt.end(), Vinf * sin(phi[i] + PI / 2.0 - AoA) + 0.5 * aeroBlade::aeroGamma + sumV);
        sumV = 0;
    }

    std::ofstream CP_outstream("output/misc/Cp.dat");
    for(int i = 0; i < numPanels; i++)
    { 
        int j = i + 1;
        if( j >= numPanels)
        {
            j = 0;
        }
        aeroCp.insert(aeroCp.end(), 1.0 - pow( (Vt[i] / (double)Vinf) , 2.0 ));
        aeroCN.insert(aeroCN.end(), -aeroCp[i] * Sj[j] * sin(phi[i] + PI / 2.0 - AoA) );
        aeroCA.insert(aeroCA.end(), -aeroCp[i] * Sj[j] * cos(phi[i] + PI / 2.0 - AoA) );

        CP_outstream << pointX[i] << " " << 1.0 - pow( (Vt[i] / (double)Vinf) , 2.0 ) << "\n";
    }
    CP_outstream.close();       

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