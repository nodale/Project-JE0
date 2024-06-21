#include <fstream>
#include <cmath>
#include <vector>

class blockMeshGen
{

std::ofstream out;
std::vector<std::vector<double>> vertices;

public:

//everything must be executed in order
void init(char* target);

void collectVertices(double x, double y, double r);

void generateVertices();

void generateEdges();

void generateBlocks();

void generateBoundaries();

private:


};