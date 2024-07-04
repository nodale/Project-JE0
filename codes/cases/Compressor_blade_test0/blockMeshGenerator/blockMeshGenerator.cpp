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
    vertices.insert(vertices.end(),{ x, y, r});
}

void blockMeshGen::generateVertices()
{
    std::string begin = "vertices\n(\n";
    out << begin;

    for(int i = 0; i < boundary.size();)
    {
        out << "( " << boundary[i][0] << " " << boundary[i][1] << " " << boundary[i][2] << " )\n";
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

    for(int i = 1; i < boundary.size();)
    {
        out << "arc " << edgeNum << " " << edgeNum + 1 << " ( " << boundary[i][0] << " " << boundary[i][1] << " " << boundary[i][2] << " )\n";
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
    for(int i = 1; i < boundary.size(); i++)
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
    out << "solid blade\n";
    
    generateFacesVerticesS(resolution);

    // for(int d = 0; d < 4; d++)
    //         {       
    //                 if(d == 0)
    //                 {
    //                 tempValue.clear();
    //                 tempValue.push_back(vertices[vertices.size() - 1]);
    //                 tempValue.push_back(vertices[vertices.size() - resolution]); 
    //                 tempValue.push_back(vertices[vertices.size() - resolution - 1]);        
    //                 }
    //                 if(d == 1)
    //                 {
    //                 tempValue.clear();
    //                 tempValue.push_back(vertices[0]);
    //                 tempValue.push_back(vertices[resolution - 1]); 
    //                 tempValue.push_back(vertices[resolution]);  
    //                 }
    //                 if(d == 2)
    //                 {
    //                 tempValue.clear();
    //                 tempValue.push_back(vertices[0]);
    //                 tempValue.push_back(vertices[resolution]); 
    //                 tempValue.push_back(vertices[resolution + 1]);  
    //                 }
    //                 if(d == 3)
    //                 {
    //                 tempValue.clear();
    //                 tempValue.push_back(vertices[vertices.size() - 1]);
    //                 tempValue.push_back(vertices[vertices.size() - resolution - 1]); 
    //                 tempValue.push_back(vertices[vertices.size() - resolution - 2]);        
    //                 }
    //                 tempValue.clear();

    //                 for(int e = 0; e < 2; e++)
    //                 {
    //                     for(int f = 0; f < 3; f++)
    //                     {
    //                         if(std::isnan(tempValue[e][f]) == 1)
    //                         {
    //                             tempValue[e][f] = 0.0;
    //                         }

    //                     }
    //                 }
                    
    //                 for(int i = 0; i < 3; i++)
    //                 {
    //                     tempLine[0][i] = tempValue[1][i] - tempValue[0][i]; 
    //                     tempLine[1][i] = tempValue[2][i] - tempValue[0][i]; 
    //                 }
                    
    //                 tempPerpendicular[0] = tempLine[0][1] * tempLine[1][2] - tempLine[0][2] * tempLine[1][1];
    //                 tempPerpendicular[1] = tempLine[0][2] * tempLine[1][0] - tempLine[0][0] * tempLine[1][2];
    //                 tempPerpendicular[2] = tempLine[0][1] * tempLine[1][0] - tempLine[0][0] * tempLine[1][1];

    //                 tempLength = sqrt( tempPerpendicular[0]*tempPerpendicular[0] + tempPerpendicular[1]*tempPerpendicular[1] + tempPerpendicular[2]*tempPerpendicular[2] );
                
    //                 for(int i = 0; i < 3; i++)
    //                 {
    //                     tempNormal[i] = fabs(tempPerpendicular[i] / tempLength);
    //                 }

    //                 out << "facet normal " << tempNormal[0] << " " << tempNormal[1] << " " << tempNormal[2] << "\n";
    //                 //out << "facet normal " << 0 << " " << 0 << " " << 0 << "\n";
    //                 out << "outer loop\n";
                    
    //                 for(int i = 0; i < 3; i++)
    //                 {
    //                 out << "vertex " << tempValue[i][0] << " " << tempValue[i][1] << " " << tempValue[i][2] << "\n";
    //                 }

    //                 out << "endloop\n" << "endfacet\n";         
                
    //         }
    std::cout << "blade is generated\n";

}

void blockMeshGen::collectBoundary(double x, double y, double r)
{
    boundary.insert(boundary.end(), { x, y, r });
}

void blockMeshGen::generateBoundary(int resolution)
{
    out.close();
    out.open("boundary.stl");

    out << "solid boundary\n";

    generateFacesVerticesB(resolution);    

    for(int d = 0; d < 2; d++)
    {       
            if(d == 9)
            {
            tempValue.clear();
            tempValue.push_back(boundary[boundary.size() - 1]);
            tempValue.push_back(boundary[boundary.size() - resolution]); 
            tempValue.push_back(boundary[boundary.size() - resolution - 1]);        
            }
            if(d == 9)
            {
            tempValue.clear();
            tempValue.push_back(boundary[0]);
            tempValue.push_back(boundary[resolution - 1]); 
            tempValue.push_back(boundary[resolution ]);  
            }
            if(d == 0)
            {
            tempValue.clear();
            tempValue.push_back(boundary[0]);
            tempValue.push_back(boundary[resolution ]); 
            tempValue.push_back(boundary[resolution + 1]);  
            }
            if(d == 1)
            {
            tempValue.clear();
            tempValue.push_back(boundary[boundary.size() - 1]);
            tempValue.push_back(boundary[boundary.size() - resolution - 1]); 
            tempValue.push_back(boundary[boundary.size() - resolution - 2]);        
            }
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
        
    }
    std::cout << "boundary is generated\n";
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

void blockMeshGen::generateFacesVerticesS(int resolution)
{
    for(int d = 0; d < vertices.size() - (resolution + 2); d++)
    {

        for(int l = 0; l < 2; l++)
        {

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
            //out << "facet normal " << 0 << " " << 0 << " " << 0 << "\n";vv
            out << "outer loop\n";
            
            for(int i = 0; i < 3; i++)
            {
            out << "vertex " << tempValue[i][0] << " " << tempValue[i][1] << " " << tempValue[i][2] << "\n";
            }

            out << "endloop\n" << "endfacet\n";
         
        
        }
    
    }
}

void blockMeshGen::generateFacesVerticesB(int resolution)
{
    bool skip = 1;
    
    for(int d = 0; d < boundary.size() - (resolution + 2); d++)
    {

        
        for(int l = 0; l < 2; l++)
        {   
            skip = 1;
            if( remainder(d + 1 + resolution / 2, ( resolution ) ) == 0 )
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[d+resolution+1]);
                skip = 0;
            }

            if( remainder(d + 2 + resolution / 2, ( resolution ) ) == 0)
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
                skip = 0;
            }
            
            if(remainder(d + 1, ( resolution ) ) == 0 )
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[d+resolution+1]);

                skip = 0;
            }

            if(remainder(d + 2, ( resolution ) ) == 0 ) 
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
                skip = 0;
            }

            if(l == 1 && skip == 1)
            {   
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
            }

            if(l == 0 && skip == 1)
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[d+resolution+1]);
            }
            

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
        }
    
    }

}

void blockMeshGen::generateInlet(int resolution)
{
    out.close();
    out.open("inlet.stl");

    out << "solid inlet\n";


        for(int d = 0; d < boundary.size() - (resolution + 3); d++)
        {     
            
            if(remainder(d + 1, ( resolution ) ) == 0 )
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
            }

            if(remainder(d + 2, ( resolution ) ) == 0 ) 
            {

                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[d+resolution+1]);
            }

            if( remainder(d + 2, ( resolution ) ) != 0 && remainder(d + 1, ( resolution ) ) != 0 )
            {
                tempValue.clear();
                continue;
            }

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
        }

        for(int l = 0; l < 2; l++)
        {   
            if(l == 0)
            {
            tempValue.clear();
            tempValue.push_back(boundary[boundary.size() - 1]);
            tempValue.push_back(boundary[boundary.size() - resolution]); 
            tempValue.push_back(boundary[boundary.size() - resolution - 1]);        
            }
            if(l == 1)
            {
            tempValue.clear();
            tempValue.push_back(boundary[0]);
            tempValue.push_back(boundary[resolution - 1]); 
            tempValue.push_back(boundary[resolution ]);  
            }

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
        }

    std::cout << "inlet is generated\n";

}

void blockMeshGen::generateOutlet(int resolution)
{
    out.close();
    out.open("outlet.stl");

    out << "solid outlet\n";

    for(int d = 0; d < boundary.size() - (resolution + 3); d++)
        {     
            
            if(remainder(d + 1 + resolution / 2, ( resolution ) ) == 0 )
            {
                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
            }

            if(remainder(d + 2 + resolution / 2, ( resolution ) ) == 0 ) 
            {

                tempValue.clear();
                tempValue.push_back(boundary[d+resolution+2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[d+resolution+1]);
            }

            if( remainder(d + 2 + resolution / 2, ( resolution ) ) != 0 && remainder(d + 1 + resolution / 2, ( resolution ) ) != 0 )
            {
                tempValue.clear();
                continue;
            }

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
        }

    std::cout << "outlet is generated\n";

}

void blockMeshGen::generateBot(int resolution)
{
    out.close();
    out.open("bot.stl");

    out << "solid bot\n";

    /* bool skip = 1;
    for(int d = 1; d <= resolution - 2; d++)
    {
        for(int l = 0; l < 2; l++)
        {
            skip = 1;
            if(d == resolution / 2 - 1)
            {
                if(l == 0)
                {
                tempValue.clear();
                tempValue.push_back(vertices[resolution / 2]);
                tempValue.push_back(boundary[resolution / 2 - 1]);
                tempValue.push_back(vertices[resolution / 2 - 1]);
                skip = 0;
                }
                if(l == 1)
                {
                tempValue.clear();
                tempValue.push_back(boundary[resolution / 2 - 1]);
                tempValue.push_back(vertices[resolution / 2]);
                tempValue.push_back(boundary[resolution / 2 - 0]); 
                skip = 0;
                }
            }
            if(l == 0 && skip == 1)
            {   
                tempValue.clear();
                tempValue.push_back(vertices[d]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
            }

            if(l == 1 && skip == 1) 
            {
                tempValue.clear();
                tempValue.push_back(vertices[d+1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(vertices[d]); 
            }

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
        }    
    } 

    for(int l = 0; l < 3; l++)
    {
        if(l == 0)
        {   
            tempValue.clear();
            tempValue.push_back(vertices[resolution - 0]);
            tempValue.push_back(boundary[resolution - 1]);
            tempValue.push_back(vertices[resolution - 1]);
        }

        if(l == 1)  
        {

            tempValue.clear();
            tempValue.push_back(vertices[1]);  
            tempValue.push_back(boundary[resolution - 1]);
            tempValue.push_back(vertices[resolution - 0]); 
        }

        if(l == 2) 
        {
            tempValue.clear();
            tempValue.push_back(boundary[resolution - 1]); 
            tempValue.push_back(vertices[1]);
            tempValue.push_back(boundary[1]);
        }

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
    }  
    std::cout << "bottom boundary is generated\n"; */

    for(int d = 0; d < resolution - 1; d++)
    {

        
        for(int l = 0; l < 2; l++)
        {  

            if(l == 1)
            {   
                tempValue.clear();
                tempValue.push_back(boundary[resolution - d - 1]);
                tempValue.push_back(boundary[d+1]);
                tempValue.push_back(boundary[d]);
            }

            if(l == 0)
            {
                tempValue.clear();
                tempValue.push_back(boundary[resolution - d - 2]);
                tempValue.push_back(boundary[d+1]); 
                tempValue.push_back(boundary[resolution - d - 1]);
            }
            

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
        }
    
    }
    std::cout << "bot boundary is generated\n";

}

void blockMeshGen::generateTop(int resolution)
{
    out.close();
    out.open("top.stl");

    out << "solid top\n";

    int limitB = boundary.size() - resolution - 0;  

    /* bool skip = 1;
    for(int d = 1; d <= resolution - 2; d++)
    {
        for(int l = 0; l < 2; l++)
        {
            skip = 1;
            if(d == resolution / 2 - 1)
            {
                if(l == 0)
                {
                tempValue.clear();
                tempValue.push_back(vertices[limitS + resolution / 2]);
                tempValue.push_back(boundary[limitB + resolution / 2 - 1]);
                tempValue.push_back(vertices[limitS + resolution / 2 - 1]);
                skip = 0;
                }
                if(l == 1)
                {
                tempValue.clear();
                tempValue.push_back(boundary[limitB + resolution / 2 - 1]);
                tempValue.push_back(vertices[limitS + resolution / 2]);
                tempValue.push_back(boundary[limitB + resolution / 2 - 0]); 
                skip = 0;
                }
            }
            if(l == 0 && skip == 1)
            {   
                tempValue.clear();
                tempValue.push_back(vertices[limitS + d]);
                tempValue.push_back(boundary[limitB + d+1]);
                tempValue.push_back(boundary[limitB + d]);
            }

            if(l == 1 && skip == 1) 
            {
                tempValue.clear();
                tempValue.push_back(vertices[limitS + d+1]);
                tempValue.push_back(boundary[limitB + d+1]);
                tempValue.push_back(vertices[limitS + d]); 
            }

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
        }    
    } 

    for(int l = 0; l < 3; l++)
    {
        if(l == 0)
        {   
            tempValue.clear();
            tempValue.push_back(vertices[limitS + resolution - 0]);
            tempValue.push_back(boundary[limitB + resolution - 1]);
            tempValue.push_back(vertices[limitS + resolution - 1]);
        }

        if(l == 1)  
        {

            tempValue.clear();
            tempValue.push_back(vertices[limitS + 1]);  
            tempValue.push_back(boundary[limitB + resolution - 1]);
            tempValue.push_back(vertices[limitS + resolution - 0]); 
        }

        if(l == 2) 
        {
            tempValue.clear();
            tempValue.push_back(boundary[limitB + resolution - 1]); 
            tempValue.push_back(vertices[limitS + 1]);
            tempValue.push_back(boundary[limitB + 1]);
        }

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
    }   */

    for(int d = 0; d < resolution - 1; d++)
    {

        
        for(int l = 0; l < 2; l++)
        {  

            if(l == 1)
            {   
                tempValue.clear();
                tempValue.push_back(boundary[limitB + resolution - d - 1]);
                tempValue.push_back(boundary[limitB + d+1]);
                tempValue.push_back(boundary[limitB + d]);
            }

            if(l == 0)
            {
                tempValue.clear();
                tempValue.push_back(boundary[limitB + resolution - d - 2]);
                tempValue.push_back(boundary[limitB + d+1]); 
                tempValue.push_back(boundary[limitB + resolution - d - 1]);
            }
            

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
        }
    
    }

    std::cout << "top boundary is generated\n";

}