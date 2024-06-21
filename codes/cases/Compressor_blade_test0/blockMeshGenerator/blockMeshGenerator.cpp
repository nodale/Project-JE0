#include "blockMeshGenerator.h"

void blockMeshGen::init(char* target)
{
    out.open(target);

    char entry[]= "FoamlFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n\n";
}

void blockMeshGen::collectVertices(double x, double y, double r)
{
    std::vector<double> temp{ x, y, r };
    vertices.push_back(temp);
}

void blockMeshGen::generateVertices()
{
    std::string begin = "vertices\n(\n";
    out << begin;

    for(int i = 0; i < vertices.size();)
    {
        out << "( " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << " )\n";
        i += 2;
        //+=2 because we want to use the vertices inbetween for arc interpolation
    }

    std::string end = ");\n\n";
    out << end;
}

void blockMeshGen::generateEdges()
{
    std::string begin = "edges\n(\n";
    out << begin;
    int edgeNum = 0;

    for(int i = 1; i < vertices.size();)
    {
        out << "arc " << edgeNum << " " << edgeNum + 1 << " ( " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << " )\n";
        edgeNum += 1;
        i += 2;
        //+=2 because we want to use the vertices inbetween for arc interpolation
    }

    std::string end = ");\n\n";
    out << end;
}

void blockMeshGen::generateBlocks()
{
    std::string begin = "blocks\n(\n";
    out << begin;
    out << "hex ( ";
    for(int i = 1; i < vertices.size(); i++)
    {
        out << i << " ";
    }

    out << " ) blade \n";
    //numbers of cells in each direction, x y r
    out << "( 20 15 10 )\n";
    //cell expansion ratios
    out << "simpleGrading ( 0.5 1 1 )\n";

    std::string end = "     );\n\n";
    out << end;
}

void blockMeshGen::generateBoundaries()
{
    std::string begin = "boundary\n(\n";
    out << begin;
    out << "    ( \n";
    out << "        blade";
    out << "        {";
    out << "            type wall";
    out << "        }";
    

    std::string end = "     );\n\n";
    out << end;
}