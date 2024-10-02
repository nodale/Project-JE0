#include "simBlade.h"


double getBladeDist(int i, int j, double radius)
{
    double dist = infoBlade::chord[i][j] / infoBlade::solidity[i][j];

    return 0.5 * dist * ( radius / infoBlade::mean_radi[i][j] );
}

void readVariable(int i, double& variable, std::string baseCondition, std::string conditionalVariable, sqlite3* db, sqlite3_stmt* stmt)
{
    int rc = 0;

    std::string baseCommand = "SELECT " + baseCondition + conditionalVariable +  " FROM designParam WHERE STAGE = ? ;";

    rc = sqlite3_prepare_v2(db, baseCommand.c_str(), -1, &stmt, nullptr);
    if (rc != SQLITE_OK) 
    {
        std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
    }
    rc = sqlite3_bind_int(stmt, 1, i + 1);
    if (rc != SQLITE_OK) 
    {
        std::cerr << "Failed to bind statement: " << sqlite3_errmsg(db) << std::endl;
    }
    rc = sqlite3_step(stmt);
    if (rc == SQLITE_ROW) 
    {
        variable = sqlite3_column_double(stmt, 0);
    }
    else if (rc == SQLITE_DONE) 
    {
        std::cerr << "No rows found for STAGE = " << i << std::endl;
    } 
    else 
    {
        std::cerr << "Error executing query: " << sqlite3_errmsg(db) << std::endl;
    }
}

void readConfig(double& disX, double& disY, double& backFat, int i, int j)
{
    sqlite3* db;
    sqlite3_stmt* stmt;
    sqlite3_open("output/database/db.db", &db);

    std::string conditionalVariable;

    int rc;

    if(j == 0)
    {
        conditionalVariable = "_rotor"; 
    }
    if(j == 1)
    {
        conditionalVariable = "_stator"; 
    }

    readVariable(i, disX, "disX", conditionalVariable, db, stmt);
    readVariable(i, disY, "disY", conditionalVariable, db, stmt);
    readVariable(i, backFat, "backFat", conditionalVariable, db, stmt);
}

void simBlade::storeInDatabaseRecursive()
{
    sqlite3* db;
    sqlite3_open("output/databse/db.db", &db);
    using namespace infoBlade;

    for(int i = 0; i < totalSize; i++)
    {
    storeInDesignDatabase(db, "tip_radi1", tip_radi[i][0], i + 1);
    storeInDesignDatabase(db, "tip_radi2", tip_radi[i][1], i + 1);

    storeInDesignDatabase(db, "hub_radi1", hub_radi[i][0], i + 1);
    storeInDesignDatabase(db, "hub_radi2", hub_radi[i][1], i + 1);
    }

    sqlite3_close(db);
}

void simBlade::generateAerofoilModel(int i, int j)
{
    char filename[] = "output/model/blade.stl"; //for blockMesh : blockMeshDict
    blockMeshGen::init(filename);

    double r, shape, extraR;
    r = 1.0;
    shape = 1.0;
    
    //defining the complex equation
    std::complex<double> z, seta, thetaC, displac, F, aerofoil, dummyF;
    double potentialC, streamC, Gamma;

    double multiplier = 1.0;
    double unmultiplier = 1.0;

    //enter the displacement of the circle

    double psiA, psiC, disX, disY, backFat;
    double tempW;

    double U = 1.0;

    using namespace infoBlade;

    double theta;
    double dr = ( tip_radi[i][j] - hub_radi[i][j] ) / resolution;
    double radius;
    double dummyScale = chord[i][j];

    readConfig(disX, disY, backFat, i, j);
    displac = std::complex( disX, disY );
    extraR = sqrt( pow( r - fabs( displac.real() ) , 2 ) + pow( displac.imag() , 2 ) );
    for(int R = 0; R < resolution; R++)
    {
        radius = hub_radi[i][j] + dr * R;
        aeroBlade::getAoA(i, j, R);
        using namespace aeroBlade;
        {
            if( (disX < 0 && disY >= 0) or (disX > 0 && disY <= 0) )
            {
                for(double gh = 1.0 * ( PI ); gh >= -1.0 * ( PI );) 
                {
                    dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                    aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);

                    blockMeshGen::collectVertices( aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0, radius - hub_radi[i][j] );

                    gh -= 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                }
            }

            if( (disX >= 0 && disY >= 0) or (disX < 0 && disY < 0))
            {
                for(double gh = -2.0 * ( PI ); gh <= 0.0 * ( PI );) 
                {
                    dummyF = std::complex( extraR * cos( gh ) / backFat + displac.real(), extraR * sin(gh) / backFat + displac.imag());
                    aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R]) / RadToDegree);
                    
                    blockMeshGen::collectVertices( aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0, radius - hub_radi[i][j] );
                    gh += 2.0 * ( PI ) / (double)infoBlade::resolution / multiplier;
                }
            }
        }

        if(R == 0 or R == resolution - 1)
        {
            if((disX >= 0 && disY >= 0) or (disX < 0 && disY < 0))
            {
                std::complex<double> dummyG;
                dummyF = std::complex( extraR * cos( -2.0 * PI ) / backFat + displac.real(), extraR * sin(-2.0 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );

                dummyF = std::complex( extraR * cos( -1.0 * PI ) / backFat + displac.real(), extraR * sin(-1.0 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                
                dummyF = std::complex( extraR * cos( -2.0 * PI ) / backFat + displac.real(), extraR * sin(-2.0 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
        

                dummyF = std::complex( extraR * cos( -1.5 * PI ) / backFat + displac.real(), extraR * sin(-1.5 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
            }
            if( (disX < 0 && disY >= 0) or (disX > 0 && disY <= 0) )
            {
                std::complex<double> dummyG;
                dummyF = std::complex( extraR * cos( PI ) / backFat + displac.real(), extraR * sin(PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );

                dummyF = std::complex( extraR * cos( 0.0 ) / backFat + displac.real(), extraR * sin(0.0) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                
                dummyF = std::complex( extraR * cos( PI ) / backFat + displac.real(), extraR * sin(PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectBrick( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
        

                dummyF = std::complex( extraR * cos( 0.5 * PI ) / backFat + displac.real(), extraR * sin(0.5 * PI) / backFat + displac.imag());
                aerofoil = joukowskyTransform(dummyF, shape, std::complex<double>(0.0, rotateAngle[i][j][R] ) / RadToDegree);
                
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 - getBladeDist(i, j, radius), radius - hub_radi[i][j] );
                blockMeshGen::collectInterpolate( 1.25 * aerofoil.real() * dummyScale / 4.0 , aerofoil.imag() * dummyScale / 4.0 + getBladeDist(i, j, radius), radius - hub_radi[i][j] );
            }
        }

        continue;

    }

    blockMeshGen::generateStl(1.0 * infoBlade::resolution);
    // blockMeshGen::generateBoundary(resolution/unmultiplier);
    // blockMeshGen::generateInlet(resolution/unmultiplier);
    // blockMeshGen::generateOutlet(resolution/unmultiplier);
    // blockMeshGen::generateBot(resolution/unmultiplier);
    // blockMeshGen::generateTop(resolution/unmultiplier);

    // blockMeshGen::getSeparationVector(getBladeDist(i, j, mean_radi[i][j]));

    blockMeshGen::generateVertices();
    blockMeshGen::generateEdges();
    blockMeshGen::generateBoundaries();
    blockMeshGen::generateBlocks();

    //blockMeshGen::generateBoundaryFile();
    blockMeshGen::generateSnappy();
    //blockMeshGen::generateSurfaceFeature();
    blockMeshGen::generateCreatePatch();
    blockMeshGen::generateControlDict();

    blockMeshGen::generateObj(resolution);
}

// int main()
// {
//     infoBlade::dataBaseSetUp();
//     infoBlade::initConditionSetUp();    
//     thermoBlade::init();
//     aeroBlade::genBlade(0, 0);
//     aeroBlade::drawDeHallersNumber();
//     aeroBlade::storeInDatabaseRecursive();
//     simBlade::storeInDatabaseRecursive();
//     for(int i = 0; i < infoBlade::totalSize; i++)
//     {
//     aeroBlade::drawEverthing(i, 0, 85);
//     }
//     simBlade::generateAerofoilModel(0, 0);

//     return 0;
// }