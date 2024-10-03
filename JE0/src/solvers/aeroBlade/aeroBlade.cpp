#include "aeroBlade.h"
#include <cmath>

//defining all variables

dVec<double> aeroBlade::I, aeroBlade::J, aeroBlade::K, aeroBlade::L;
sVec<double> aeroBlade::Vn, aeroBlade::Vt, aeroBlade::lambda, aeroBlade::phi, aeroBlade::Sj;
sVec<double> aeroBlade::pointX;
sVec<double> aeroBlade::pointY;
sVec<double> aeroBlade::midPointX;
sVec<double> aeroBlade::midPointY;
double aeroBlade::aeroGamma, aeroBlade::aeroCl;
double A, B, Cn, Dn, Ct, Dt, E;

//reads the aerofoilConfig.dat file
void readConfig(double& disX, double& disY, double& backFat, int j)
{
    std::string temp1;

    std::ifstream input("input/aerofoilConfig.dat");

    std::getline(input, temp1);
    disX = std::stod(temp1);
    std::getline(input, temp1);

    if(j == 1)
    {
        disY = -std::stod(temp1);
    }
    if(j == 0)
    {
        disY = std::stod(temp1);
    }

    std::getline(input, temp1);
    backFat = std::stod(temp1);
}

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

    //the following is for printing the sum of lambda and the gamma
    // double tempSum = 0;
    // for(int i = 0; i < C.size(); i++)
    // {
    //     tempSum += C[i];
    // }
    // std::cout << tempSum << std::endl;
    // std::cout << aeroBlade::aeroGamma << std::endl;
    return C;
}

std::complex<double> aeroBlade::joukowskyTransform(std::complex<double> z, double shape, std::complex<double> thetaC)
{
    return exp(thetaC) * ( z + pow(shape,2) / z) ;  
}

void storeAerofoilConfig(int i, int j, double disX, double disY, double backFat)
{
    sqlite3* db;
    sqlite3_open("output/database/db.db", &db);
    std::string text, rotorOrStator;

    if(j == 0)
    {
        rotorOrStator = "_rotor";
    }
    if(j == 1)
    {
        rotorOrStator = "_stator";
    }

    text = "disX" + rotorOrStator;
    infoBlade::storeInDesignDatabase(db, text, disX, i + 1);

    text = "disY" + rotorOrStator;
    infoBlade::storeInDesignDatabase(db, text, disY, i + 1);

    text = "backFat" + rotorOrStator;
    infoBlade::storeInDesignDatabase(db, text, backFat, i + 1);

    sqlite3_close(db);
}

void aeroBlade::genBlade(int i, int j)
{
    std::complex<double> dummyF;
    std::complex<double> displac, aerofoil;

    pointX.clear();
    pointY.clear();

    //default values for the aerofoil configuration
    double r, shape, extraR;
    r = 1.0;
    shape = 1.0;
    double multiplier = 1.0;
    double psiA, psiC, disX, disY, backFat;
    double U = 1.0;

    readConfig(disX, disY, backFat, j);
    //storeAerofoilConfig(i, j, disX, disY, backFat);

    std::ofstream output("output/misc/shape.dat");
    
    displac = std::complex<double>(disX, disY);
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );
    using namespace aeroBlade;
    {
        if( (disX < 0 && disY >= 0) or (disX > 0 && disY <= 0) )
        {
            int i = 0;
            for(double gh = 1.0 * ( PI ); gh >= -1.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, 0.0) / RadToDegree);
                
                pointX.insert(pointX.begin() + i, 0.25 * ( aerofoil.real() + 2.0 ));
                pointY.insert(pointY.begin() + i, 0.25 * aerofoil.imag());           

                // pointX[i] = 0.25 * ( aerofoil.real() + 2.0 );
                // pointY[i] = 0.25 * aerofoil.imag();     

                gh -= 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }
        }

        if( (disX >= 0 && disY >= 0) or (disX < 0 && disY < 0))
        {
            int i = 0;
            
            for(double gh = -2.0 * ( PI ); gh <= 0.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, 0.0) / RadToDegree);
                
                pointX.insert(pointX.begin() + i, 0.25 * ( aerofoil.real() + 2.0 ));
                pointY.insert(pointY.begin() + i, 0.25 * aerofoil.imag());              

                gh += 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }
        }

        //std::cout << dist << " " << smallestPoint << std::endl;
        //for viewing the aerofoil shape, stored in output/misc/shape.dat
        double pivotPoint;
        if( disX < 0)
        {
            pivotPoint = std::min_element(pointX.begin(), pointX.end())[0];
            auto index = std::find(pointX.begin(),pointX.end(), pivotPoint);
            size_t dist = std::distance(pointX.begin(), index);
            
            for(int i = (int)dist; i < pointX.size(); i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
            for(int i = 0; i < (int)dist; i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
        }

        if( disX >= 0)
        {
            pivotPoint = std::max_element(pointX.begin(), pointX.end())[0];
            auto index = std::find(pointX.begin(),pointX.end(), pivotPoint);
            size_t dist = std::distance(pointX.begin(), index);     

            for(int i = (int)dist; i < pointX.size(); i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
            for(int i = 0; i < (int)dist; i++)
            {
                output << pointX[i] << " " << pointY[i] << std::endl;
            }
        }
        output.close();

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
        input.close();

        for(int R = 0; R < infoBlade::resolution; R++)
        {
            getAoA(i, j, R);
        }
    }
    
}

//draws the blade for stage [I, J], for a specific R, put it in a for loop
void aeroBlade::drawEverthing(int I, int j, int R)
{
    std::complex<double> dummyF;
    std::complex<double> displac, aerofoil;

    //default values for the aerofoil configuration
    double r, shape, extraR;
    r = 1.0;
    shape = 1.0;
    double multiplier = 1.0;
    double psiA, psiC, disX, disY, backFat;
    double U = 1.0;
    psiA = 0.022;
    psiC = 0.012;
    backFat = 1.02;
    disX = psiA;// - psiC * cos( 0 );
    disY = 0.034;

    readConfig(disX, disY, backFat, j);

    aeroBlade::pointX.clear();
    aeroBlade::pointY.clear();

    aeroBlade::pointX.resize(infoBlade::resolution);
    aeroBlade::pointY.resize(infoBlade::resolution);

    displac = std::complex<double>(disX, disY);
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );
    using namespace aeroBlade;
    {
        if( (disX < 0 && disY >= 0) or (disX > 0 && disY <= 0) )
        {
            int i = 0;
            for(double gh = 1.0 * ( PI ); gh >= -1.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, infoBlade::incidenceAngle[I][j][R]) / RadToDegree);
                
                pointX.insert(pointX.begin() + i, 0.25 * ( aerofoil.real() + 2.0 ));
                pointY.insert(pointY.begin() + i, 0.25 * aerofoil.imag());           

                gh -= 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }
        }

        if( (disX >= 0 && disY >= 0) or (disX < 0 && disY < 0))
        {
            int i = 0;
            
            for(double gh = -2.0 * ( PI ); gh <= 0.0 * ( PI );) 
            {
                dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, infoBlade::incidenceAngle[I][j][R]) / RadToDegree);
                
                pointX.insert(pointX.begin() + i, 0.25 * ( aerofoil.real() + 2.0 ));
                pointY.insert(pointY.begin() + i, 0.25 * aerofoil.imag());              

                gh += 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                i++;
            }
        }

        //std::cout << dist << " " << smallestPoint << std::endl;
        //for viewing the aerofoil shape, stored in output/misc/shape.dat

        std::ofstream output("output/misc/shapeAll.dat", std::ios_base::app);     
        for(int i = 0; i < infoBlade::resolution; i++)
        {
            output << pointX[i] + I * 4.0 + j * 2.0 << " " << pointY[i] << std::endl;
        }

        output << "\n";

        pointX.clear();
        pointY.clear();

        
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

void calculateRequiredCl()
{
    using namespace infoBlade;
    for(int i = 0; i < totalSize; i++)
    {
    double theta, yValue;

    for(int r = 0; r <= resolution; r++)
    {
    
    // if(j == 0)
    // {
    // theta = cos (beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    // }
    // if(j == 1)
    // {
    // theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    // }

    //get loss coefficients, not extremely accurate but it's alright
    lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i][1] / cos( alpha[i+1][0][r] / RadToDegree );

    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    //pressure loss is in Kpascal
    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

    //batchAnalysis[i] << theta << " " <<  liftCoefficient[i][j][r] << "\n";

    //efficiency analysis
    //double dummyEfficiency = 1 - ( lossCoefficient[i][0][r] / pow ( cos( alpha[i+1][0][r] / RadToDegree ), 2 ) + lossCoefficient[i][1][r] / pow ( cos( beta[i][0][r] / RadToDegree ), 2 ) ) * pow( phi[i] , 2 ) / ( 2 * psi[i] );
    //batchAnalysis[i] << r << " " <<  dummyEfficiency << "\n";
    }
    }
}

void aeroBlade::getAoA(int i, int j, int R)
{
    double dischargeY1, dischargeY2, dischargeX1, dischargeX2;
    using namespace infoBlade;

    dischargeAngle[i][j][R] = atan2(  aeroBlade::pointY[resolution - 2] - aeroBlade::pointY[resolution - 3] ,  aeroBlade::pointX[1] - aeroBlade::pointX[2] ) * RadToDegree;

    if( j == 0)
    {
    rotateAngle[i][j][R] = beta[i][1][R] - dischargeAngle[i][j][R];
    incidenceAngle[i][j][R] = beta[i][0][R] - rotateAngle[i][j][R];
    }

    if( j == 1)
    {
    rotateAngle[i][j][R] = alpha[i+1][0][R] - dischargeAngle[i][j][R];  
    incidenceAngle[i][j][R] = alpha[i][1][R] - rotateAngle[i][j][R];
    }

    //std::cout << "incidence angle : " << incidenceAngle[i][j][R] << std::endl;
}

void aeroBlade::storeInDatabaseRecursive()
{
    sqlite3* db;
    sqlite3_open("output/databse/db.db", &db);
    using namespace infoBlade;

    for(int i = 0; i < totalSize; i++)
    {
    storeInAeroDatabase(db, "meanAlpha1", meanAlpha[i][0], i + 1);
    storeInAeroDatabase(db, "meanAlpha2", meanAlpha[i][1], i + 1);
    storeInAeroDatabase(db, "meanBeta1", meanBeta[i][0], i + 1);
    storeInAeroDatabase(db, "meanBeta2", meanBeta[i][1], i + 1);
    storeInAeroDatabase(db, "PressureRatio", PR[i], i + 1);
    storeInAeroDatabase(db, "Area1", Area[i][0], i + 1);
    storeInAeroDatabase(db, "Area2", Area[i][1], i + 1);
    storeInAeroDatabase(db, "Area3", Area[i][2], i + 1);
    storeInAeroDatabase(db, "chord1", chord[i][0], i + 1);
    storeInAeroDatabase(db, "chord2", chord[i][1], i + 1);
    storeInAeroDatabase(db, "solidity1", solidity[i][0], i + 1);
    storeInAeroDatabase(db, "solidity2", solidity[i][1], i + 1);
    }

    sqlite3_close(db);
}

void aeroBlade::drawDeHallersNumber()
{
    std::ofstream out;
    std::string path;
    double theta;

    std::string gnuCommand = "#!/usr/bin/gnuplot -persist\nplot ";
    std::string generateExecutable;

    for(int j = 0; j < 2; j++)
    {   
        for(int i = 0; i < infoBlade::totalSize; i++)
        {
            path = "output/batchData/" + std::to_string(i) + "_" + std::to_string(j) + ".dat";
            out.open(path);
            using namespace infoBlade;
            for(int r = 0; r < infoBlade::resolution; r++)
            {
                if(j == 0)
                {
                    theta = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
                }
                if(j == 1)          
                {
                    theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
                }

                lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
                lossCoefficient[i][1][r] = 0.014 * solidity[i][j] / cos( alpha[i+1][0][r] / RadToDegree );
                
                Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
                Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

                pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
                pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

                liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
                liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

                out << theta << " " << infoBlade::liftCoefficient[i][j][r] << "\n";
            }
            out.close();

            gnuCommand += "'output/batchData/" + std::to_string(i) + "_" + std::to_string(j) + ".dat"+ "' with linespoints linetype -1 linewidth 2 lc rgb '#6194B',";
        }
    }

    generateExecutable = "echo \"" + gnuCommand + "\" > bin/drawDeHaller.sh";

    std::system("touch bin/drawDeHaller.sh\n chmod +x bin/drawDeHaller.sh");
    std::system(generateExecutable.c_str());

}

void aeroBlade::drawDeHallersNumber_v2()
{
    std::ofstream out;
    std::string path;
    double theta, deltaBeta;

    std::string gnuCommand = "#!/usr/bin/gnuplot -persist\nplot ";
    std::string generateExecutable;

    for(int j = 0; j < 2; j++)
    {   
        for(int i = 0; i < infoBlade::totalSize; i++)
        {
            path = "output/batchData/" + std::to_string(i) + "_" + std::to_string(j) + ".dat";
            out.open(path);
            using namespace infoBlade;
            for(int r = 0; r < infoBlade::resolution; r++)
            {
                if(j == 0)
                {
                    theta = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
                    deltaBeta = beta[i][0][r]/RadToDegree - beta[i][1][r]/RadToDegree;
                    liftCoefficient[i][0][r] = 2 * deltaBeta * work[i] * R[i] / chord[i][j];
                }
                if(j == 1)          
                {
                    theta = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
                    deltaBeta = beta[i][1][r]/RadToDegree - beta[i + 1][0][r]/RadToDegree;
                    liftCoefficient[i][1][r] = 2 * deltaBeta * work[i] * ( 1.0 - R[i] ) / chord[i][j];
                }

                out << theta << " " << infoBlade::liftCoefficient[i][j][r] << "\n";
            }
            out.close();

            gnuCommand += "'output/batchData/" + std::to_string(i) + "_" + std::to_string(j) + ".dat"+ "' with linespoints linetype -1 linewidth 2 lc rgb '#6194B',";
        }
    }

    generateExecutable = "echo \"" + gnuCommand + "\" > bin/drawDeHaller.sh";

    std::system("touch bin/drawDeHaller.sh\n chmod +x bin/drawDeHaller.sh");
    std::system(generateExecutable.c_str());

}

int aeroBlade::findCombinationAlpha(int sampleSize, int maxTries)
{
    std::uniform_int_distribution<> distr2(20, 40); //omega1
    std::uniform_int_distribution<> distr3(-70, 70); //alpha1 angles
    float dir1, dir2;
    double theta1, theta2, smallestGradient;
    std::random_device rd1;
    dVec<int> suitableAlphas;
    dVec<int> successfulAlphas;
    sVec<double> accumulatedSmallestGradient;
    sVec<double> localSmallestGradient;
    sVec<double> stageSmallestGradient;
    dVec<double> tempGradient;
    sVec<double> suitableOmega1;
    int tempOmega1;
    int j = 0;
    int tries = 0;
    using namespace infoBlade;

    tempGradient.resize(totalSize);
    suitableAlphas.resize(sampleSize);

    bool firstAttempt = true;

    dir1 = 1.8;
    dir2 = -1.8;
    
    for(int nklf = 0; nklf < sampleSize; nklf++) //sample size
    {

    for(int i = 0; i < totalSize ; i++)
    {

    for(int r = 0; r < resolution; r++)    
    {
    
    tries += 1;
    if(tries >= maxTries)
    {
        std::cout << "no design fits the given RPM, increase the maxTries or change the RPM\n";
        std::cout << "number of combinations found so far : " << nklf << std::endl;
        return 0;
    }
    
    while(true)
    {               
   
    theta1 = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
    theta2 = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
    
    lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
    lossCoefficient[i][1][r] = 0.014 * solidity[i][j] / cos( alpha[i+1][0][r] / RadToDegree );
    
    Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
    Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

    pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
    pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

    liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
    liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

    std::cout << solidity[i][0] << " " << solidity[i][1] << " " << theta1 << " " << theta2 << " " << liftCoefficient[i][0][r] << " " << liftCoefficient[i][1][r] << std::endl;

    if
    (
    theta1 >= 0.73 && theta2 >= 0.73 && liftCoefficient[i][0][r] <= dir1 && liftCoefficient[i][0][r] >= dir2 && solidity[i][0] <= 5.0
     && firstAttempt == false && liftCoefficient[i][1][r] <= dir1 && liftCoefficient[i][1][r] >= dir2 && solidity[i][1] <= 5.0
    )
    {
        if(r > 0)
        {
        tempGradient[i].insert(tempGradient[i].end(), liftCoefficient[i][j][r] - liftCoefficient[i][j][r-1]);
        }
        if(r == resolution -1 && i == totalSize - 1)
        {
            for( int g = 0; g < totalSize; g++)
            {
                stageSmallestGradient.insert(stageSmallestGradient.end(), std::max_element(tempGradient[g].begin(), tempGradient[g].end())[0]);
                suitableAlphas[nklf].insert(suitableAlphas[nklf].end(), meanAlpha[g][0]);
                suitableOmega1.insert(suitableOmega1.end(), tempOmega1);
                tempGradient[g].clear();
                //std::cout << stageSmallestGradient[j] << std::endl;
            }
            localSmallestGradient.insert(localSmallestGradient.end(), std::max_element(stageSmallestGradient.begin(), stageSmallestGradient.end())[0]);
            
            //std::cout << "combination no " << nklf << " passes\n";      

            stageSmallestGradient.clear();

            firstAttempt = true;
            tries = 0;
        }
        break;
    }
    else                                        
    {          
        firstAttempt = false; 
        tempGradient[i].clear();     
        suitableAlphas[nklf].clear();

        //int tempA = distr3(rd1) * 2;        
        for(int m = 0; m < lowSize; m++)
        {
            //PR[m] = tempB;
            //omega1 = distr2(rd1) * 200;
            tempOmega1 = omega1;
            meanAlpha[m][0] = distr3(rd1);
        }
        //tempA = distr3(rd1) * 2;
        for(int m = lowSize; m < totalSize; m++)
        {
            //PR[m] = pow( 14 / pow( tempB, 3 ) , 1.0 / 8.0 );
            //somega1 = distr2(rd1) * 200;
            meanAlpha[m][0] = distr3(rd1);
        }
        
        //std::cout << "trying the following combination: " << omega1 << " //" << omega2 << " ";

        thermoBlade::calculateThermoVariables();

        //std::cout << "\n";

        r = 0;
        i = 0;
    }

    }

    }
    }

    //accumulatedSmallestGradient.insert(accumulatedSmallestGradient.end(), std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0]);

    }
    // for(int l = 0; l < sampleSize; l++)
    // {
    //     std::cout << localSmallestGradient[l] << std::endl;
    // }
    smallestGradient = std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0];
    auto index = std::find(localSmallestGradient.begin(),localSmallestGradient.end(), smallestGradient);
    size_t dist = std::distance(localSmallestGradient.begin(), index);

    std::cout << "The following combination has been stored : ";

    for(int c = 0; c < totalSize; c++)
    {
        meanAlpha[c][0] = suitableAlphas[dist][c];
        std::cout << suitableAlphas[dist][c] << " ";
    }

    thermoBlade::calculateThermoVariables();
    storeInDatabaseRecursive();

    std::cout << "\n";
    //std::cout << omega1 << std::endl;
    return 0;
}

int aeroBlade::findCombinationFull(int sampleSize, int maxTries)
{
    std::uniform_int_distribution<> distr2(2, 700); //omega1
    std::uniform_int_distribution<> distr3(-70, 70); //alpha1 angles
    float dir1, dir2;
    double theta1, theta2, smallestGradient;
    std::random_device rd1;
    dVec<int> suitableAlphas;
    dVec<int> successfulAlphas;
    sVec<double> accumulatedSmallestGradient;
    sVec<double> localSmallestGradient;
    sVec<double> stageSmallestGradient;
    dVec<double> tempGradient;
    sVec<double> suitableOmega1;
    sVec<double> suitableOmega2;
    int tempOmega1 = infoBlade::omega1;
    int tempOmega2 = infoBlade::omega2;
    int j = 0;
    int tries = 0;
    using namespace infoBlade;

    tempGradient.resize(totalSize);
    suitableAlphas.resize(sampleSize);

    bool firstAttempt = true;

    dir1 = 1.8;
    dir2 = -1.8;
    
    for(int nklf = 0; nklf < sampleSize; nklf++) //sample size
    {

        for(int i = 0; i < totalSize ; i++)
        {

            for(int r = 0; r < resolution; r++)    
            {
            
                tries += 1;
                if(tries >= maxTries)
                {
                    std::cout << "no design fits the given RPM, increase the maxTries or change the RPM\n";
                    std::cout << "number of combinations found so far : " << nklf << std::endl;
                    return 0;
                }
                
                while(true)
                {               
                theta1 = cos(beta[i][0][r]/RadToDegree) / cos(beta[i][1][r]/RadToDegree); 
                theta2 = cos(alpha[i+1][0][r]/RadToDegree) / cos(alpha[i][1][r]/RadToDegree); 
                
                lossCoefficient[i][0][r] = 0.014 * solidity[i][0] / cos( beta[i][1][r] / RadToDegree );
                lossCoefficient[i][1][r] = 0.014 * solidity[i][j] / cos( alpha[i+1][0][r] / RadToDegree );
                
                Mach[i][0][r] = VX[0] / ( cos( beta[i][1][r] / RadToDegree ) * pow( Temperature[i][1] * 287 * gamma , 0.5 ) );
                Mach[i][1][r] = VX[0] / ( cos( alpha[i+1][0][r] / RadToDegree ) * pow( Temperature[i][2] * 287 * gamma , 0.5 ) );

                pressureLoss[i][0][r] = 0.000005 * rho[i][1] * lossCoefficient[i][0][r] * pow( VX[0] / cos( beta[i][1][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][0][r] , 2 ) , gamma / ( gamma - 1 ) );
                pressureLoss[i][1][r] = 0.000005 * rho[i][2] * lossCoefficient[i][1][r] * pow( VX[0] / cos( alpha[i+1][0][r] / RadToDegree ) , 2 ) * pow( 1 + 0.5 * ( gamma - 1 ) * pow( Mach[i][1][r] , 2 ) , gamma / ( gamma - 1 ) );

                liftCoefficient[i][0][r] = 2.0 / solidity[i][0] * ( tan( beta[i][0][r] / RadToDegree ) - tan( beta[i][1][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][0][r] * sin( 0.5 * ( tan( beta[i][0][r] / RadToDegree ) + tan( beta[i][1][r] / RadToDegree ) ) ) / ( rho[i][1] * pow( 0.5 * ( VX[0] / cos( beta[i][0][r] / RadToDegree ) + VX[0] / cos( beta[i][1][r] / RadToDegree ) ) , 2 ) * solidity[i][0] );
                liftCoefficient[i][1][r] = 2.0 / solidity[i][1] * ( tan( alpha[i][1][r] / RadToDegree ) - tan( alpha[i+1][0][r] / RadToDegree ) ) * cos( atan( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) ) - 2 * pressureLoss[i][1][r] * sin( 0.5 * ( tan( alpha[i][1][r] / RadToDegree ) + tan( alpha[i+1][0][r] / RadToDegree ) ) ) / ( rho[i][2] * pow( 0.5 * ( VX[1] / cos( alpha[i][1][r] / RadToDegree ) + VX[1] / cos( alpha[i+1][0][r] / RadToDegree ) ) , 2 ) * solidity[i][1] );

                //std::cout << solidity[i][0] << " " << solidity[i][1] << " " << theta << " " << liftCoefficient[i][0][r] << " " << liftCoefficient[i][1][r] << std::endl;

                if
                (
                theta1 >= 0.73 && theta2 >= 0.73 && liftCoefficient[i][0][r] <= dir1 && liftCoefficient[i][0][r] >= dir2 && solidity[i][0] <= 5.0
                && firstAttempt == false && liftCoefficient[i][1][r] <= dir1 && liftCoefficient[i][1][r] >= dir2 && solidity[i][1] <= 5.0
                )
                {
                    if(r > 0)
                    {
                    tempGradient[i].insert(tempGradient[i].end(), liftCoefficient[i][j][r] - liftCoefficient[i][j][r-1]);
                    }
                    if(r == resolution - 1 && i == totalSize - 1)
                    {
                        for( int g = 0; g < totalSize; g++)
                        {
                            stageSmallestGradient.insert(stageSmallestGradient.end(), std::max_element(tempGradient[g].begin(), tempGradient[g].end())[0]);
                            suitableAlphas[nklf].insert(suitableAlphas[nklf].end(), meanAlpha[g][0]);
                            suitableOmega1.insert(suitableOmega1.end(), tempOmega1);
                            suitableOmega2.insert(suitableOmega1.end(), tempOmega2);
                            tempGradient[g].clear();
                            //std::cout << stageSmallestGradient[j] << std::endl;
                        }
                        localSmallestGradient.insert(localSmallestGradient.end(), std::max_element(stageSmallestGradient.begin(), stageSmallestGradient.end())[0]);
                        
                        //std::cout << "combination no " << nklf << " passes\n";      

                        stageSmallestGradient.clear();

                        firstAttempt = true;
                        tries = 0;
                    }
                    break;
                }
                else                                        
                {          
                    firstAttempt = false; 
                    tempGradient[i].clear();     
                    suitableAlphas[nklf].clear();

                    //int tempA = distr3(rd1) * 2;        
                    for(int m = 0; m < lowSize; m++)
                    {
                        //PR[m] = tempB;
                        omega1 = distr2(rd1) * 20;
                        tempOmega1 = omega1;
                        meanAlpha[m][0] = distr3(rd1);
                    }
                    //tempA = distr3(rd1) * 2;
                    for(int m = lowSize; m < totalSize; m++)
                    {
                        //PR[m] = pow( 14 / pow( tempB, 3 ) , 1.0 / 8.0 );
                        omega2 = distr2(rd1) * 20;
                        tempOmega2 = omega2;
                        meanAlpha[m][0] = distr3(rd1);
                    }
                    
                    //std::cout << "trying the following combination: " << omega1 << " //" << omega2 << " ";

                    thermoBlade::calculateThermoVariables();

                    //std::cout << "\n";

                    r = 0;
                    i = 0;
                }

                }

            }
        }

    //accumulatedSmallestGradient.insert(accumulatedSmallestGradient.end(), std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0]);

    }
    // for(int l = 0; l < sampleSize; l++)
    // {
    //     std::cout << localSmallestGradient[l] << std::endl;
    // }
    smallestGradient = std::min_element(localSmallestGradient.begin(), localSmallestGradient.end())[0];
    auto index = std::find(localSmallestGradient.begin(),localSmallestGradient.end(), smallestGradient);
    size_t dist = std::distance(localSmallestGradient.begin(), index);

    std::cout << "The following combination has been stored : ";

    for(int c = 0; c < totalSize; c++)
    {
        meanAlpha[c][0] = suitableAlphas[dist][c];
        std::cout << suitableAlphas[dist][c] << " ";
    }

    omega1 = suitableOmega1[dist];

    thermoBlade::calculateThermoVariables();
    storeInDatabaseRecursive();

    std::cout << "with omega1 : " << omega1 << " and omega2 : " << omega2 << std::endl;
    std::cout << "\n";
    return 0;
}

void aeroBlade::sourceVortexPanelMethodAngle(int i, int j, int r)
{
    using namespace aeroBlade;
    double temp;
    dVec<double> avg;
    dVec<double> matrixA, inverseMatrixA;
    sVec<double> matrixB, aeroCp, aeroCN, aeroCA, aeroCm, aeroCL;

    double AoA = -infoBlade::incidenceAngle[i][j][r];
    double Vinf = 1.0; //infoBlade::VX[0];

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

    //std::ofstream CL_outstream("output/misc/Cl.dat");
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
        aeroCL.insert(aeroCL.end(), aeroCN[i] * cos(AoA) - aeroCA[i] * sin(AoA));
        //CL_outstream << pointX[i] << " " << aeroCL[i] << "\n";
        CP_outstream << pointX[i] << " " << aeroCp[i] << "\n";
    }
    //CL_outstream.close();    
    CP_outstream.close();    

    double sumLength = 0.0;
    for(int i = 0; i < numBoundaries; i++)
    {
        sumLength += Sj[i];         
    }   
    aeroBlade::aeroCl = sumLength * 2 * aeroBlade::aeroGamma;

    std::cout << "total lift coefficient : " << aeroBlade::aeroCl << std::endl;

}   

void aeroBlade::sourceVortexPanelMethodNoAngle(int i, int j, int r)
{
    using namespace aeroBlade;
    double temp;
    dVec<double> avg;
    dVec<double> matrixA, inverseMatrixA;
    sVec<double> matrixB, aeroCp, aeroCN, aeroCA, aeroCm, aeroCL;

    double AoA = 0.0;
    double Vinf = 1.0; //infoBlade::VX[0];

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

    //std::ofstream CL_outstream("output/misc/Cl.dat");
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
        aeroCL.insert(aeroCL.end(), aeroCN[i] * cos(AoA) - aeroCA[i] * sin(AoA));
        //CL_outstream << pointX[i] << " " << aeroCL[i] << "\n";
        CP_outstream << pointX[i] << " " << aeroCp[i] << "\n";
    }
    //CL_outstream.close();    
    CP_outstream.close();    

    double sumLength = 0.0;
    for(int i = 0; i < numBoundaries; i++)
    {
        sumLength += Sj[i];         
    }   
    aeroBlade::aeroCl = sumLength * 2 * aeroBlade::aeroGamma;

    std::cout << "total lift coefficient : " << aeroBlade::aeroCl << std::endl;

}

//only for testing  
// int main()
// {       
//     infoBlade::dataBaseSetUp();
//     infoBlade::initConditionSetUp();    
//     thermoBlade::init();
//     aeroBlade::drawDeHallersNumber();
//     storeInDatabaseRecursive();

//     int r = infoBlade::resolution / 2;

//     //aeroBlade::findCombinationAlpha(5,10000000);

//     for(int j = 0; j < 2; j++)
//     {
//         for(int i = 0; i < infoBlade::totalSize; i++)
//         {
//         aeroBlade::genBlade(i, j);
//         getAoA(i, j, r);
//         aeroBlade::sourceVortexPanelMethod(i, j ,r); 
//         drawEverthing(i, j, r);                          
//         }
//     }
//     calculateRequiredCl();   

//     return 0;
// }