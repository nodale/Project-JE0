#include "../infoBlade/infoBlade.h"

const char* argv = 
{
"CREATE TABLE thermoBlade ("
"   Temperature1     REAL,"
"   Temperature2     REAL,"
"   Temperature3     REAL,"
"   TemperatureStag1 REAL,"
"   TemperatureStag2 REAL,"
"   TemperatureStag3 REAL,"
"   Pressure1        REAL,"
"   Pressure2        REAL,"
"   Pressure3        REAL,"
"   PressureStag1    REAL,"
"   PressureStag2    REAL,"
"   PressureStag3    REAL,"
"   rho1             REAL,"
"   rho2             REAL,"
"   rho3             REAL,"
"   psi              REAL,"
"   phi              REAL,"
"   a                REAL,"
"   b                REAL,"
"   Wr               REAL,"
"   Ws               REAL,"
"   efficiency1      REAL,"
"   efficiency2      REAL,"  
"   efficiency3      REAL,"  
");"

"CREATE TABLE aeroBlade ("
"   Mach0            REAL,"
"   Mach1            REAL,"
"   meanAlpha1       REAL,"
"   meanAlpha2       REAL,"
"   meanBeta1        REAL,"
"   meanBeta2        REAL,"
"   PressureRatio    REAL,"
"   Area1            REAL,"
"   Area2            REAL,"
"   Area3            REAL,"
"   chord1           REAL,"
"   chord2           REAL,"
"   solidity1        REAL,"
"   solidity2        REAL,"
");"

"CREATE TABLE designParam ("
"   tip_radi1        REAL,"
"   tip_radi2        REAL,"
"   tip_radi3        REAL,"
"   hub_radi1        REAL,"
"   hub_radi2        REAL,"
"   hub_radi3        REAL,"
"   numBlades1       REAL,"
"   numBlades2       REAL,"
");"
};

//this does NOTHING
static int callback(void* data, int argc, char** argv, char** ColName)
{
    return 0;
}

void checkRC(int rc)
{
    if( rc )
    {
        std::cout << "ERROR : can't open database\n";
    }
}

//seets up the tables for all parameters
bool infoBlade::dataBaseSetUp()
{
    sqlite3* db;
    int rc;
    char* ErrorMsg = 0;

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    rc = sqlite3_exec(db, argv, callback, 0, &ErrorMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "ERROR : failed setting up dataBase\n";
        return 1;
    }
    else
    {
        return 0;
    }

    sqlite3_close(db);
}

//reads initial conditions (T1, P1, and rho0) and put them into dataBase
bool infoBlade::initConditionSetUp()
{
    sqlite3* db;
    int rc;
    char* ErrorMsg = 0;

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    

    sqlite3_close(db);
}