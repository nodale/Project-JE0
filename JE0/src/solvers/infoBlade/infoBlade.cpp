#include "infoBlade.h"
#include <cstdio>

//TODO
//embed the other solvers code with SQLite 

//static variables definitions used in this source file
int infoBlade::totalSize;     
int infoBlade::lowSize;
int infoBlade::highSize;
double infoBlade::T1;
double infoBlade::P1;
double infoBlade::rho1;    

//other static variables
dVec<double> hub_radi;
dVec<double> mean_radi;
double omega1;
double omega2;              
dVec<double> v;
sVec<double> work;
dVec<double> diffusion;
dVec<double> numBlades;
int lowSize;
int highSize;
int totalSize;
dVec<std::vector<double>> Mach;
dVec<std::vector<double>> alpha;
dVec<std::vector<double>> beta;
dVec<double> meanAlpha;
dVec<double> meanBeta;
double Mach0;
double resolution;
sVec<double> PR;
sVec<double> R;
dVec<double> Area;
dVec<double> chord;
dVec<double> maxCam;       
dVec<double> maxCamPos;
dVec<std::vector<double>> liftCoefficient;
dVec<std::vector<double>> rotateAngle;
dVec<std::vector<double>> incidenceAngle;
dVec<std::vector<double>> AoA;
double dummyMaxCam, dummyMaxCamPos, dummyChord, dummyLiftCoefficient;
dVec<std::vector<double>> lossCoefficient;
dVec<std::vector<double>> pressureLoss;
dVec<std::vector<double>> dischargeAngle;
dVec<double> solidity;
dVec<double> Temperature;
dVec<double> TemperatureStag;
dVec<double> Pressure;
dVec<double> PressureStag;
dVec<double> rho;
sVec<double> psi;
sVec<double> phi;
sVec<double> a;
sVec<double> b;
sVec<double> Wr;
sVec<double> Ws;
dVec<double> efficiency;


const char* argv = 
{
"CREATE TABLE thermoBlade ("
"   STAGE            INT     PRIMARY KEY,"
"   Temperature1     REAL    NOT NULL DEFAULT 0.0,"

"   Temperature2     REAL    NOT NULL DEFAULT 0.0,"
"   Temperature3     REAL    NOT NULL DEFAULT 0.0,"
"   TemperatureStag1 REAL    NOT NULL DEFAULT 0.0,"
"   TemperatureStag2 REAL    NOT NULL DEFAULT 0.0,"
"   TemperatureStag3 REAL    NOT NULL DEFAULT 0.0,"
"   Pressure1        REAL    NOT NULL DEFAULT 0.0,"
"   Pressure2        REAL    NOT NULL DEFAULT 0.0,"
"   Pressure3        REAL    NOT NULL DEFAULT 0.0,"
"   PressureStag1    REAL    NOT NULL DEFAULT 0.0,"
"   PressureStag2    REAL    NOT NULL DEFAULT 0.0,"
"   PressureStag3    REAL    NOT NULL DEFAULT 0.0,"
"   rho1             REAL    NOT NULL DEFAULT 0.0,"
"   rho2             REAL    NOT NULL DEFAULT 0.0,"
"   rho3             REAL    NOT NULL DEFAULT 0.0,"
"   psi              REAL    NOT NULL DEFAULT 0.0,"
"   phi              REAL    NOT NULL DEFAULT 0.0,"
"   a                REAL    NOT NULL DEFAULT 0.0,"
"   b                REAL    NOT NULL DEFAULT 0.0,"
"   Wr               REAL    NOT NULL DEFAULT 0.0,"
"   Ws               REAL    NOT NULL DEFAULT 0.0,"
"   efficiency1      REAL    NOT NULL DEFAULT 0.0,"
"   efficiency2      REAL    NOT NULL DEFAULT 0.0,"  
"   efficiency3      REAL    NOT NULL DEFAULT 0.0"  
");"

"CREATE TABLE aeroBlade ("
"   STAGE            INT     PRIMARY KEY,"
"   Mach0            REAL    NOT NULL DEFAULT 0.0,"
"   Mach1            REAL    NOT NULL DEFAULT 0.0,"
"   meanAlpha1       REAL    NOT NULL DEFAULT 0.0,"
"   meanAlpha2       REAL    NOT NULL DEFAULT 0.0,"
"   meanBeta1        REAL    NOT NULL DEFAULT 0.0,"
"   meanBeta2        REAL    NOT NULL DEFAULT 0.0,"
"   PressureRatio    REAL    NOT NULL DEFAULT 0.0,"
"   Area1            REAL    NOT NULL DEFAULT 0.0,"
"   Area2            REAL    NOT NULL DEFAULT 0.0,"
"   Area3            REAL    NOT NULL DEFAULT 0.0,"
"   chord1           REAL    NOT NULL DEFAULT 0.0,"
"   chord2           REAL    NOT NULL DEFAULT 0.0,"
"   solidity1        REAL    NOT NULL DEFAULT 0.0,"
"   solidity2        REAL    NOT NULL DEFAULT 0.0"
");"

"CREATE TABLE designParam ("
"   STAGE            INT     PRIMARY KEY,"
"   tip_radi1        REAL    NOT NULL DEFAULT 0.0,"
"   tip_radi2        REAL    NOT NULL DEFAULT 0.0,"
"   tip_radi3        REAL    NOT NULL DEFAULT 0.0,"
"   hub_radi1        REAL    NOT NULL DEFAULT 0.0,"
"   hub_radi2        REAL    NOT NULL DEFAULT 0.0,"
"   hub_radi3        REAL    NOT NULL DEFAULT 0.0,"
"   numBlades1       REAL    NOT NULL DEFAULT 0.0,"
"   numBlades2       REAL    NOT NULL DEFAULT 0.0"
");"
};

const char* argc =
{
"INSERT INTO thermoBlade (Temperature1, Pressure1, rho1) VALUES ();"
};

//this does NOTHING
static int callback(void* data, int argc, char** argv, char** ColName)
{
    return 0;
}

static int insertStages(void* data, int argc, char** argv, char** ColName)
{
    std::cout << "A\n";
    return 0;
}

void checkRC(int rc)
{
    if( rc )
    {
        std::cout << "ERROR : can't open database\n";
    }
}

void openCurrentInput(std::ifstream& input)
{
    std::filesystem::path newPath;
    std::filesystem::path currentPath = getCurrentPath();
    newPath = currentPath / "input" / "initCondInput.dat";
    input.open(newPath);
}

//seets up the tables for all parameters
bool infoBlade::dataBaseSetUp()
{
    std::string tempL, tempH;
    std::ifstream input;
    openCurrentInput(input);
    sqlite3* db;
    int rc;
    char* ErrorMsg = 0;
    sqlite3_stmt* stmt;
    const char* rowInput[3] = 
    {
        "INSERT INTO thermoBlade (STAGE) VALUES (?);",
        "INSERT INTO aeroBlade (STAGE) VALUES (?);",
        "INSERT INTO designParam (STAGE) VALUES (?);",

    };

    std::getline(input, tempL);
    std::getline(input, tempH);
    infoBlade::totalSize = std::stoi(tempL) + std::stoi(tempH); 

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    rc = sqlite3_exec(db, argv, callback, NULL, &ErrorMsg);
    if (rc != SQLITE_OK)            
    {
        std::cout << "ERROR : failed setting up dataBase ; " << ErrorMsg << std::endl;
    }

    for(int j = 0; j < 3; j++)
    {
        rc = sqlite3_prepare_v2(db, rowInput[j], -1, &stmt, NULL);
        if(rc != SQLITE_OK)
        {
            std::cout << "ERROR : statement preparation failed ; " << sqlite3_errmsg(db) << std::endl;
        }
        for(int i = 1; i <= totalSize; i++)
        {
            sqlite3_bind_int(stmt, 1, i);
            rc = sqlite3_step(stmt);
            if(rc != SQLITE_DONE)
            {
                std::cout << "ERROR : SQLite step failed ; " << rc << std::endl;
            }
            sqlite3_reset(stmt);
        }
    }
    sqlite3_finalize(stmt);

    sqlite3_close(db);

    return 0;
}

//reads initial conditions (T1, P1, and rho0) and put them into dataBase
bool infoBlade::initConditionSetUp()
{
    std::string temp1;
    sqlite3* db;
    int rc;
    char* ErrorMsg = 0;

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    std::ifstream input;
    openCurrentInput(input);
    std::getline(input, temp1);
    infoBlade::lowSize = std::stoi(temp1);
    std::getline(input, temp1);
    infoBlade::highSize = std::stoi(temp1);
    std::getline(input, temp1);
    infoBlade::T1 = std::stod(temp1);
    std::getline(input, temp1);
    infoBlade::P1 = std::stod(temp1);
    std::getline(input, temp1);
    infoBlade::rho1 = std::stod(temp1);    

    sqlite3_close(db);

    return 0;
}               

//only for testing
int main()
{
    infoBlade test;

    test.dataBaseSetUp();
    test.initConditionSetUp();

    std::cout << infoBlade::rho1 << std::endl;

    return 0;
}