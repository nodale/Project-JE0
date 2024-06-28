#include "blockMeshGenerator.h"
#include <cmath>

void blockMeshGen::init(char* target)
{
    out.open(target);

    //enable the following two for blockMesh
    //char entry[]= "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n}\n\n\n";
    //out << entry;
}

void blockMeshGen::collectVertices(double x, double y, double r)
{
    std::vector<double> temp{ x, y, r };
    vertices.insert(vertices.end(),temp);
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
    out << "        blade\n";
    out << "        {\n";
    out << "            type wall\n";
    out << "        }\n";
    out << "    )\n";
    

    std::string end = "     );\n\n";
    out << end;
}

void blockMeshGen::generateStl(int resolution)
{
    //calculate normal
    std::vector<std::vector<double>> tempValue;
    double tempLine[2][3];
    double tempPerpendicular[3];
    double tempNormal[3];
    double tempLength;
    int quotient; 

    int dummyRes = resolution + 1;
    out << "solid blade\n";
    std::vector<std::vector<double>> tempVertices;

    // tempVertices.push_back(vertices[resolution + 2]);
    // tempVertices.push_back(vertices[resolution]);
    // tempVertices.push_back(vertices[resolution+1]);
    // tempVertices.push_back(vertices[1]);
    

    for(int d = 0; d < vertices.size() - (resolution + 2); d++)
    {
        for(int l = 0; l < 2; l++)
        {
            // if( remainder( d - resolution, dummyRes) != 0 )
            // {
            
            // if( remainder( d - resolution, dummyRes) == 0 )
            // {
            //     if(l == 0)
            //     {
            //         tempValue.clear();
            //         tempValue.push_back(tempVertices[3]);
            //         tempValue.push_back(tempVertices[0]);
            //         tempValue.push_back(tempVertices[2]);
            //     }

            //     if(l == 1)
            //     {
            //         tempValue.clear();
            //         tempValue.push_back(tempVertices[1]);
            //         tempValue.push_back(tempVertices[0]);
            //         tempValue.push_back(tempVertices[3]);
            //     }
            // }

            // if( remainder( d - 1, resolution) != 0  )
            // {
                if(l == 0)
                {   
                    tempValue.clear();
                    tempValue.push_back(vertices[d+resolution+1]);
                    tempValue.push_back(vertices[d+1]);
                    tempValue.push_back(vertices[d]);
                }

                if(l == 1)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[d+resolution+ 2]);
                    tempValue.push_back(vertices[d+1]); 
                    tempValue.push_back(vertices[d+resolution+1]);
                }
            //} 
            


            //need to fix this
            for(int e = 0; e < 2; e++)
            {
                for(int f = 0; f < 3; f++)
                {
                    if(std::isnan(tempValue[e][f]) == 1)
                    {
                        tempValue[e][f] = 0.0;
                    }

                }
            }
            
            for(int i = 0; i < 3; i++)
            {
                tempLine[0][i] = tempValue[1][i] - tempValue[0][i]; 
                tempLine[1][i] = tempValue[2][i] - tempValue[0][i]; 
            }
            
            tempPerpendicular[0] = tempLine[0][1] * tempLine[1][2] - tempLine[0][2] * tempLine[1][1];
            tempPerpendicular[1] = tempLine[0][2] * tempLine[1][0] - tempLine[0][0] * tempLine[1][2];
            tempPerpendicular[2] = tempLine[0][1] * tempLine[1][0] - tempLine[0][0] * tempLine[1][1];

            tempLength = sqrt( tempPerpendicular[0]*tempPerpendicular[0] + tempPerpendicular[1]*tempPerpendicular[1] + tempPerpendicular[2]*tempPerpendicular[2] );
        
            for(int i = 0; i < 3; i++)
            {
                tempNormal[i] = fabs(tempPerpendicular[i] / tempLength);
            }

            out << "facet normal " << tempNormal[0] << " " << tempNormal[1] << " " << tempNormal[2] << "\n";
            //out << "facet normal " << 0 << " " << 0 << " " << 0 << "\n";
            out << "outer loop\n";
            
            for(int i = 0; i < 3; i++)
            {
            out << "vertex " << tempValue[i][0] << " " << tempValue[i][1] << " " << tempValue[i][2] << "\n";
            }

            out << "endloop\n" << "endfacet\n";

            

            // if( remainder( d, dummyRes) == 0 )
            // {
            //     quotient = d / dummyRes;
    
            //     tempVertices.clear();
            //     tempVertices.push_back(vertices[d + resolution + 1]);
            //     tempVertices.push_back(vertices[d]);
            //     tempVertices.push_back(vertices[d + resolution + 2]);
            //     tempVertices.push_back(vertices[d+1]);    
            // }            
        
        }
    
    }

vertices.clear();    
}

void blockMeshGen::clear()
{
    vertices.clear();
}

void blockMeshGen::generateSnappy()
{
    std::ofstream output("CPD/snappyHexMeshDict");

    output << "Foamfile\n";
    output << "{\n";
    output << " version     2.0;\n";
    output << " format      ascii;\n";
    output << " class       dictionary;\n";
    output << " object      autoHexMeshDict;\n";
    output << "}\n";

    output << "castellatedMesh  true;\n"; 

    output << "geometry\n";
    output << "{\n";
    output << " blade.stl\n";
    output << " {\n";
    output << "     type triSurfaceMesh;\n";
    output << "     name blade;\n";
    output << " }\n";
    output << "}\n";   
}