#include "../infoBlade/infoBlade.h"

//TODO
//find out why the stages are not inserted

const char* argv = 
{
"CREATE TABLE thermoBlade ("
"   STAGE            INT     PRIMARY KEY,"
"   Temperature1     REAL    NOT NULL UNIQUE,"
"   Temperature2     REAL    NOT NULL UNIQUE,"
"   Temperature3     REAL    NOT NULL UNIQUE,"
"   TemperatureStag1 REAL    NOT NULL UNIQUE,"
"   TemperatureStag2 REAL    NOT NULL UNIQUE,"
"   TemperatureStag3 REAL    NOT NULL UNIQUE,"
"   Pressure1        REAL    NOT NULL UNIQUE,"
"   Pressure2        REAL    NOT NULL UNIQUE,"
"   Pressure3        REAL    NOT NULL UNIQUE,"
"   PressureStag1    REAL    NOT NULL UNIQUE,"
"   PressureStag2    REAL    NOT NULL UNIQUE,"
"   PressureStag3    REAL    NOT NULL UNIQUE,"
"   rho1             REAL    NOT NULL UNIQUE,"
"   rho2             REAL    NOT NULL UNIQUE,"
"   rho3             REAL    NOT NULL UNIQUE,"
"   psi              REAL    NOT NULL UNIQUE,"
"   phi              REAL    NOT NULL UNIQUE,"
"   a                REAL    NOT NULL UNIQUE,"
"   b                REAL    NOT NULL UNIQUE,"
"   Wr               REAL    NOT NULL UNIQUE,"
"   Ws               REAL    NOT NULL UNIQUE,"
"   efficiency1      REAL    NOT NULL UNIQUE,"
"   efficiency2      REAL    NOT NULL UNIQUE,"  
"   efficiency3      REAL    NOT NULL UNIQUE"  
");"

"CREATE TABLE aeroBlade ("
"   STAGE            INT     PRIMARY KEY,"
"   Mach0            REAL    NOT NULL UNIQUE,"
"   Mach1            REAL    NOT NULL UNIQUE,"
"   meanAlpha1       REAL    NOT NULL UNIQUE,"
"   meanAlpha2       REAL    NOT NULL UNIQUE,"
"   meanBeta1        REAL    NOT NULL UNIQUE,"
"   meanBeta2        REAL    NOT NULL UNIQUE,"
"   PressureRatio    REAL    NOT NULL UNIQUE,"
"   Area1            REAL    NOT NULL UNIQUE,"
"   Area2            REAL    NOT NULL UNIQUE,"
"   Area3            REAL    NOT NULL UNIQUE,"
"   chord1           REAL    NOT NULL UNIQUE,"
"   chord2           REAL    NOT NULL UNIQUE,"
"   solidity1        REAL    NOT NULL UNIQUE,"
"   solidity2        REAL    NOT NULL UNIQUE"
");"

"CREATE TABLE designParam ("
"   STAGE            INT     PRIMARY KEY,"
"   tip_radi1        REAL    NOT NULL UNIQUE,"
"   tip_radi2        REAL    NOT NULL UNIQUE,"
"   tip_radi3        REAL    NOT NULL UNIQUE,"
"   hub_radi1        REAL    NOT NULL UNIQUE,"
"   hub_radi2        REAL    NOT NULL UNIQUE,"
"   hub_radi3        REAL    NOT NULL UNIQUE,"
"   numBlades1       REAL    NOT NULL UNIQUE,"
"   numBlades2       REAL    NOT NULL UNIQUE"
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

    sqlite3_stmt* stmt;
    const char* rowInput[3] = 
    {
        "INSERT INTO thermoBlade (STAGE) VALUES (?);",
        "INSERT INTO aeroBlade (STAGE) VALUES (?);",
        "INSERT INTO designParam (STAGE) VALUES (?);",

    };

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    rc = sqlite3_exec(db, argv, callback, NULL, &ErrorMsg);
    if (rc != SQLITE_OK)            
    {
        std::cout << "ERROR : failed setting up dataBase ; " << ErrorMsg << std::endl;
        return 1;
    }
    else
    {
        return 0;
    }

    for(int j =0; j < 3; j++)
    {
        sqlite3_reset(stmt);
        rc = sqlite3_prepare16_v2(db, rowInput[j], -1, &stmt, NULL);
        if(rc != SQLITE_OK)
        {
            std::cout << "ERROR : statement preparation failed\n";
        }
        for(int i = 0; i < totalSize; i++)
        {
            sqlite3_bind_int(stmt, 1, i);
            rc = sqlite3_step(stmt);
            if(rc != SQLITE_DONE)
            {
                std::cout << "ERROR : SQLite step failed\n";
            }
        }
    }
    sqlite3_finalize(stmt);

    sqlite3_close(db);

    return 0;
}

//reads initial conditions (T1, P1, and rho0) and put them into dataBase
bool infoBlade::initConditionSetUp()
{
    sqlite3* db;
    int rc;
    char* ErrorMsg = 0;

    rc = sqlite3_open("output/database/db.db", &db);
    checkRC(rc);

    std::filesystem::path currentPath = getCurrentPath();
    std::ifstream input(currentPath);

    

    sqlite3_close(db);

    return 0;
}

//only for testing
int main()
{
    infoBlade test;

    test.dataBaseSetUp();

    return 0;
}