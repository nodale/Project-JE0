#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>

class blockMeshGen
{

std::ofstream out;
std::vector<std::vector<double>> vertices;
std::vector<std::vector<double>> boundary;

std::vector<std::vector<double>> tempValue;
double tempLine[2][3];
double tempPerpendicular[3];
double tempNormal[3];
double tempLength;

public:

void init(char* target);
void collectVertices(double x, double y, double r);
void collectBoundary(double x, double y, double r);
//everything must be executed in order, only for blockMesh
void generateVertices();
void generateEdges();
void generateBlocks();
void generateBoundaries();

//generate an stl file
void generateStl(int resolution);
void generateBoundary(int resolution);
void generateInlet(int resolution);
void generateOutlet(int resolution);
void generateTop(int resolution);
void generateBot(int resolution);


//still in work, for snappyHexMesh
void generateSnappy();

void clear();

private:

void generateFacesVerticesS(int resolution);
void generateFacesVerticesB(int resolution);
void generateFacesVerticesBCombined(int resolution);



};