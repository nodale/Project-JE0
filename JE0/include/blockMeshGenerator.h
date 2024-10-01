#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

namespace blockMeshGen
{

extern std::ofstream out;
extern std::ofstream bmD;
extern std::ofstream bmB;
extern std::ofstream CFDFolder;

extern std::vector<std::vector<double>> vertices;
extern std::vector<std::vector<double>> boundary;
extern std::vector<std::vector<double>> brick;
extern std::vector<std::vector<double>> interpolatePoints;

extern double separationVector;

extern std::vector<std::vector<double>> tempValue;
extern double tempLine[2][3];
extern double tempPerpendicular[3];
extern double tempNormal[3];
extern double tempLength;

void init(char* target);
void collectVertices(double x, double y, double r);
void collectBoundary(double x, double y, double r);
void collectBrick(double x, double y, double r);
void collectInterpolate(double x, double y, double r);

//everything must be executed in order, only for blockMesh
void generateVertices();
void generateEdges();
void generateBlocks();
void generateBoundaries();
void generateBoundaryFile();
void generateSurfaceFeature();

//generate an stl file
void generateStl(int resolution);
void generateBoundary(int resolution);
void generateInlet(int resolution);
void generateOutlet(int resolution);
void generateTop(int resolution);
void generateBot(int resolution);

//generate an obj file
void generateObj(int resolution);

void getSeparationVector(double input);

//still in work, for snappyHexMesh
void generateSnappy();
void generateCreatePatch();
void generateControlDict();

void clear();


void generateFacesVerticesS(int resolution);
void generateFacesVerticesB(int resolution);
void generateFacesVerticesBCombined(int resolution);

void printCoordinates();
void printCoordinatesObj(std::ofstream& obj);
};