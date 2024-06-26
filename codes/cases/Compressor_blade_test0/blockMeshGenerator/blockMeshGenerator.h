#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>

class blockMeshGen
{

std::ofstream out;
std::vector<std::vector<double>> vertices;

public:

void init(char* target);
void collectVertices(double x, double y, double r);
//everything must be executed in order, only for blockMesh
void generateVertices();
void generateEdges();
void generateBlocks();
void generateBoundaries();

//generate an stl file
void generateStl(int resolution);

void clear();

private:


};