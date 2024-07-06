#include "blockMeshGenerator.h"
#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>

void blockMeshGen::init(char* target)
{
    out.open(target);
    bmD.open("systemFile/blockMeshDict");
    bmB.open("boundary/boundary");


    //enable the following two for blockMesh
    char entry[]= "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      blockMeshDict;\n}\n\n\n\n\n";
    bmD << entry;
    bmD << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
}

void blockMeshGen::collectVertices(double x, double y, double r)
{
    vertices.insert(vertices.end(),{ x, y, r});
}

void blockMeshGen::collectBrick(double x, double y, double r)
{
    brick.insert(brick.end(),{ x, y, r});
}

void blockMeshGen::collectInterpolate(double x, double y, double r)
{
    interpolatePoints.insert(interpolatePoints.end(),{ x, y, r});
}

void blockMeshGen::printCoordinates()
{
    // std::vector<std::vector<float>> tempCheck(tempValue.size());
    // std::transform(tempValue.begin(), tempValue.end(), tempCheck.begin(), [](const std::vector<double>& vec) {
    //                    std::vector<float> vecF(vec.size());
    //                    std::transform(vec.begin(), vec.end(), vecF.begin(),
    //                                   [](double val) { return static_cast<float>(val); });
    //                    return vecF;
    //                });

    // if(tempCheck[0] == tempCheck[1] or tempCheck[0] == tempCheck[2] or tempCheck[1] == tempCheck[2])
    // {
    //     std::cout << "collapsed a face\n";
    //     return;
    // }

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

void blockMeshGen::generateVertices()
{
    std::string begin = "vertices\n(\n";
    bmD << begin;

    for(int i = brick.size() - 1; i >= 0;)
    {
        bmD << "    ( " << brick[i][0] << " " << brick[i][1] << " " << brick[i][2] << " )\n";
        i -= 1;
        //+=2 because we want to use the vertices inbetween for arc interpolation
    }

    std::string end = ");\n\n";
    bmD << end;
}

void blockMeshGen::generateEdges()
{
    std::string begin = "edges\n(\n";
    bmD << begin;
    int edgeNum = 0;

    // for(int i = 1; i < boundary.size();)
    // {
    //     bmD << "arc " << edgeNum << " " << edgeNum + 1 << " ( " << boundary[i][0] << " " << boundary[i][1] << " " << boundary[i][2] << " )\n";
    //     edgeNum += 1;
    //     i += 2;
    //     //+=2 because we want to use the vertices inbetween for arc interpolation
    // }
    
    // for(int i = 0; i <= 4; i++)
    // {
    //     bmD << "arc " << 0 << " " << 1 << " ( " << interpolatePoints[i][0] << " " << interpolatePoints[i][1] << " " << interpolatePoints[i][2] << " )\n";   
    // }

    bmD << "arc " << 6 << " " << 7 << " ( " << interpolatePoints[0][0] << " " << interpolatePoints[0][1] << " " << interpolatePoints[0][2] << " )\n";
    bmD << "arc " << 4 << " " << 5 << " ( " << interpolatePoints[1][0] << " " << interpolatePoints[1][1] << " " << interpolatePoints[1][2] << " )\n";

    bmD << "arc " << 2 << " " << 3 << " ( " << interpolatePoints[2][0] << " " << interpolatePoints[2][1] << " " << interpolatePoints[2][2] << " )\n";
    bmD << "arc " << 0 << " " << 1 << " ( " << interpolatePoints[3][0] << " " << interpolatePoints[3][1] << " " << interpolatePoints[3][2] << " )\n";


    std::string end = ");\n\n";
    bmD << end;
}

void blockMeshGen::generateBlocks()
{
    std::string begin = "blocks\n(\n";
    bmD << begin;
    bmD << "hex ( 0 1 2 3 4 5 6 7";

    bmD << " ) boundary \n";
    //numbers of cells in each direction, x y r
    bmD << "( 70 70 50 )\n";
    //cell expansion ratios
    bmD << "simpleGrading ( 2 2 2 )\n";

    std::string end = ");\n\n";
    bmD << end;
}

void blockMeshGen::generateBoundaries()
{
    std::string begin = "boundary\n(\n";
    bmD << begin;
    bmD << "        inlet\n"
     "        {\n"
     "            type outletInlet;\n"
     "            faces\n"
     "            (\n"
     "                   ( 0 3 4 7 )\n"
     "            );\n"
     "         }\n\n"
     "        outlet\n"
     "        {\n"
     "            type outletInlet;\n"
     "            faces\n"
     "            (\n"
     "                   ( 1 2 5 6 )\n"
     "            );\n"
     "         }\n\n"
     "        casing\n"
     "        {\n"
     "            type wall;\n"
     "            faces\n"
     "            (\n"
     "                   ( 0 1 2 3 )\n"
     "                   ( 4 5 6 7 )\n"
     "            );\n"
     "         }\n\n"
     "        cyclic1\n"
     "        {\n"
     "            type              wall;\n"
    //  "            inGroups          1(cyclic);\n"
    //  "            matchTolerance    0.002;\n"
    //  "            neighbourPatch    cyclic2;\n"
    //  "            transform         rotational;\n"
    //  "            rotationAxis      (0 0 1);\n"
    //  "            rotationCentre    (0 0 0);\n"
     "            faces\n"
     "            (\n"
     "                   ( 0 1 4 5 )\n"
     "            );\n"
     "         }\n\n"
     "        cyclic2\n"
     "        {\n"
     "            type              wall;\n"
    //  "            inGroups          1(cyclic);\n"
    //  "            matchTolerance    0.002;\n"
    //  "            neighbourPatch    cyclic1;\n"
    //  "            transform         rotational;\n"
    //  "            rotationAxis      (0 0 1);\n"
    //  "            rotationCentre    (0 0 0);\n"
     "            faces\n"
     "            (\n"
    "                   ( 2 3 6 7 )\n"
    "            );\n"
    "         }\n\n";
    

    std::string end = ");\n\n";
    bmD << end;
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
    out.open("stl/boundary.stl");

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

void blockMeshGen::generateBoundaryFile()
{
    char entry[]= "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       polyBoundaryMesh;\n    object      blockMeshDict;\n}\n\n\n";
    bmB << entry;

    bmB << "5\n"
    "(";
    bmB << 
    "    inlet\n"
    "    {\n"
    "       type            outletIntlet;\n"
    "       outletValue     (2, 2, 2);\n"
    "       value           2;\n"
    "    }\n"
    "    outlet\n"
    "    {\n"
    "       type            outletInlet;\n"
    "       outletValue     (2,2,2);\n"
    "       value           2;\n"
    "    }\n"
    "    casing\n"
    "    {\n"
    "       type            wall;\n"
    "       nFaces          50;\n"
    "       startFace       10325;\n"
    "    }\n"
    "    cyclic1\n"
    "    {\n"
    "       type            cyclic;\n"
    "       nFaces          50;\n"
    "       startFace       10350;\n"
    "       matchTolerance  0.15;\n"
    "       neighbourPatch  cyclic2;\n"
    "    }\n"
    "    cyclic2\n"
    "    {\n"
    "       type            cyclic;\n"
    "       nFaces          50;\n"
    "       startFace       10375;\n"
    "       matchTolerance  0.15;\n"
    "       neighbourPatch  cyclic1;\n"
    "    }\n"
    ")";



}
void blockMeshGen::generateSnappy()
{
    std::ofstream output("systemFile/snappyHexMeshDict");

    /* output << "FoamFile\n"
    "{\n"
    "   version     2.0;\n"
    "   format      ascii;\n"
    "   class       dictionary;\n"
    "   object      snappyHexMeshDict;\n"
    "}\n\n"
    "\n"  
    "\n"  
    "castellatedMesh  true;\n" 
    "snap             true;\n" 
    "addLayers        false;\n\n" 
    "geometry\n"
    "{\n"
    "   blade.stl\n"
    "   {\n"
    "       type triSurfaceMesh;\n"
    "       name aerofoil;\n"
    "\n"  
    "       regions\n"  
    "       {\n"
    "           blade\n"  
    "           {\n"  
    "               name bladeSolid;\n"  
    "            }\n"  
    "       }\n"  
    "       locationInMesh (0 0 0.02);\n"  
    "   }\n"
    "}\n"
    "\n"  
    "castellatedMeshControls\n"  
    "{\n"
    "   maxGlobalCells 20000000;\n"  
    "   maxLocalCells 1000000;\n"  
    "   minRefinementCells 2;\n"  
    "   nCellsBetweenLevels 2;\n"
    "   maxLoadUnbalance 0.10;\n"  
    "\n"  
    "   features\n"  
    "   (\n"  
    "       {\n"  
    "       file \"blade.eMesh\";\n"  
    "       level 2;\n"  
    "       }\n"  
    "   );\n"  
    "\n"  
    "\n"  
    "   refinementSurfaces\n"  
    "   {\n"  
    "       aerofoil\n"  
    "       {\n"  
    "           level (1 1);\n"  
    "           regions\n"  
    "           {\n"  
    "               blade\n"
    "               {\n"  
    "                   level (1 1);\n"  
    "                   patchInfo\n"  
    "                   {\n"
    "                       type wall;\n"  
    "                   }\n"  
    "               }\n"  
    "           }\n"  
    "           faceZone    aerofoil;\n"  
    "           faceType    boundary;\n"  
    "           cellZone    aerofoil;\n"  
    "           cellZoneInside  inside;\n"  
    "\n"  
    "\n"  
    "\n"  
    "       }\n"  
    "   }\n"  
    "\n" 
    "\n"  
    "   refinementRegions\n"  
    "   {\n"   
    "   }\n"  
    "\n"  
    "   resolveFeatureAngle 10;\n"  
    "\n"    
    "   allowFreeStandingZoneFaces true;\n"
    "}\n"  
    "\n"
    "   cleanMesh true;\n"
    "\n"
    "snapControls\n"
    "{\n"
    "   nSmoothPatch 3;\n"
    "   tolerance 2.0;\n"
    "   nSolveIter 300;\n"
    "   nRelaxIter 3;\n"
    "   nFeatureSnapIter 10;\n"
    "\n"  
    "   implicitFeatureSnap     false;\n"  
    "   explicitFeatureSnap     true;\n"  
    "   multiRegionFeatureSnap  false;\n"  
    "\n"  
    "}\n"
    "\n"
    "\n"

    "\n"
    "addLayersControls\n"  
    "{\n"  
    "    layers\n"  
    "    {\n"  
    "       bladesolid\n"
    "       {\n"  
    "           nSurfaceLayers 3;\n"  
    "       }\n"  
    "    }\n"  
    "\n"  
    "\n"  
    "\n"  
    "\n"  
    "   nBufferCellsNoExtrude 0;\n"  
    "   maxFaceThicknessRatio 0.5;\n"  
    "   nGrow 0;\n"  
    "\n"  
    "   nSmoothThickness 8;\n"  
    "   nSmoothNormals 2;\n"
    "   meanMedianAxisAngle 90;\n"  
    "   nSmoothSurfaceNormals 1;\n"  
    "   maxThicknessToMedialRatio   0.3;\n"  
    "\n"  
    "\n"  
    "\n"  
    "\n"  
    "   slipFeatureAngle 30;\n"  
    "   featureAngle 130;\n"  
    "   minThickness 0.0004;\n" 
    "   finalLayerThickness 0.6;\n"  
    "   expansionRatio 1.02;\n"  
    "   relativeSizes true;\n"  
    "\n"  
    "   nLayerIter 50;\n"  
    "   nRelaxIter 5;\n"  
    "   nLayerIter 20;\n"  
    "   cleanMesh true;\n"  
    "\n"  
    "\n"  
    "}\n"  
    "\n"  
    "\n"
    "meshQualityControls\n"
    "{\n"
    "   maxConcave 160;\n"
    "   maxNonOrtho 75;\n"
    "   minVol 1e-13;\n"
    "   minDeterminant 0.001;\n"
    "   minFaceWeight 0.05;\n"
    "   minVolRatio 0.01;\n"
    "   minTetQuality 1e-06;\n"  
    "   minArea 1e-30;\n"  
    "   minVol 1e-30;\n"  
    "   maxInternalSkewness 4;\n"  
    "   maxBoundarySkewness 20;\n"  
    "   minTwist -1;\n"  
    "   minTriangleTwist -1;\n"  
    "   nSmoothScale 4;\n"  
    "   errorReduction 0.75;\n"  
    "\n"  
    "\n"  
    "}\n"
    "\n"
    "mergeTolerance 1e-02;\n"
    "debug 3;\n"
    "\n"; */

    output << "FoamFile\n"
    "{\n"
    "   version     2.0;\n"
    "   format      ascii;\n"
    "   class       dictionary;\n"
    "   object      snappyHexMeshDict;\n"
    "}\n\n"
    "\n"  
    "\n"  
    "castellatedMesh  true;\n" 
    "snap             false;\n" 
    "addLayers        false;\n\n" 
    "geometry\n"
    "{\n"
    "   blade\n"
    "   {\n"
    "       type triSurfaceMesh;\n"
    "       file \"blade.stl\";\n"
    "   }\n"
    "}\n"
    "\n"  
    "castellatedMeshControls\n"  
    "{\n"
    "   maxGlobalCells 20000000;\n"  
    "   maxLocalCells 1000000;\n"  
    "   minRefinementCells 100;\n"  
    "   nCellsBetweenLevels 2;\n"
    "   maxLoadUnbalance 0.10;\n"  
    "\n"  
    "   features\n"  
    "   (\n"   
    "   );\n"  
    "\n"  
    "\n"  
    "   refinementSurfaces\n"  
    "   {\n"  
    "       blade\n"  
    "       {\n"  
    "           level (3 8);\n"  
    "           patchInfo\n"  
    "           {\n"
    "               type wall;\n"  
    "           }\n"  
    "       }\n"  
    "   }\n"  
    "\n" 
    "\n"  
    "   refinementRegions\n"  
    "   {\n"   
    "   }\n"  
    "\n"  
    "   resolveFeatureAngle 10;\n"  
    "   locationInMesh  (0 0 0);\n"    
    "   allowFreeStandingZoneFaces true;\n"
    "}\n"  
    "\n"
    "   cleanMesh true;\n"
    "\n"
    "snapControls\n"
    "{\n"
    "   nSmoothPatch 3;\n"
    "   tolerance 2.0;\n"
    "   nSolveIter 30;\n"
    "   nRelaxIter 5;\n"
    "   nFeatureSnapIter 10;\n"
    "\n"  
    "   implicitFeatureSnap     false;\n"  
    "   explicitFeatureSnap     true;\n"  
    "   multiRegionFeatureSnap  false;\n"  
    "\n"  
    "}\n"
    "\n"
    "\n"

    "\n"
    "addLayersControls\n"  
    "{\n"  
    "    layers\n"  
    "    {\n"  
    "       bladesolid\n"
    "       {\n"  
    "           nSurfaceLayers 3;\n"  
    "       }\n"  
    "    }\n"  
    "\n"  
    "\n"  
    "\n"  
    "\n"  
    "   nBufferCellsNoExtrude 0;\n"  
    "   maxFaceThicknessRatio 0.5;\n"  
    "   nGrow 0;\n"  
    "\n"  
    "   nSmoothThickness 8;\n"  
    "   nSmoothNormals 2;\n"
    "   meanMedianAxisAngle 90;\n"  
    "   nSmoothSurfaceNormals 1;\n"  
    "   maxThicknessToMedialRatio   0.3;\n"  
    "\n"  
    "\n"  
    "\n"  
    "\n"  
    "   slipFeatureAngle 30;\n"  
    "   featureAngle 130;\n"  
    "   minThickness 0.0004;\n" 
    "   finalLayerThickness 0.6;\n"  
    "   expansionRatio 1.02;\n"  
    "   relativeSizes true;\n"  
    "\n"  
    "   nLayerIter 50;\n"  
    "   nRelaxIter 5;\n"  
    "   nLayerIter 20;\n"  
    "   cleanMesh true;\n"  
    "\n"  
    "\n"  
    "}\n"  
    "\n"  
    "\n"
    "meshQualityControls\n"
    "{\n"
    "   maxConcave          180;\n"
    "   maxNonOrtho         180;\n"
    "   minVol              1e-13;\n"
    "   minDeterminant      0.001;\n"
    "   minFaceWeight       0.05;\n"
    "   minVolRatio         0.01;\n"
    "   minVolCollapseRatio 0.1;\n"
    "   minEdgeLength       -1;\n"
    "   minTetQuality       -1e-15;\n"  
    "   minArea             -1;\n"  
    "   minVol              -1e-15;\n"  
    "   maxInternalSkewness 4;\n"  
    "   maxBoundarySkewness 20;\n"  
    "   minTwist            -1;\n"  
    "   minTriangleTwist    -1;\n"  
    "   nSmoothScale        4;\n"  
    "   errorReduction      0.75;\n"  
    "\n"  
    "\n"  
    "}\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "mergeTolerance 1e-02;\n"
    "debug 3;\n"
    "\n";


}
void blockMeshGen::generateControlDict()
{
    std::ofstream cd("systemFile/controlDict");
    cd << "FoamFile\n"
    "{\n"
    "   version     2.0;\n"
    "   format      ascii;\n"
    "   class       dictionary;\n"
    "   object      controlDict;\n"
    "}\n\n"
    "\n"
    "application        icoFoam;\n"
    "\n"
    "startFrom          startTime;\n"
    "\n"
    "startTime          0.0;\n"
    "\n"
    "stopAt             endTime;\n"
    "\n"
    "endTime            0.0;\n"
    "\n"
    "writeControl       timeStep;\n"
    "\n"
    "writeInterval      1;\n"
    "\n"
    "deltaT             0.01;\n"
    "\n"
    "purgeWrite         0;\n"
    "\n"
    "writeFormat        ascii;\n"
    "\n"
    "writePrecision     6.0;\n"
    "\n"
    "writeCompression   off;\n"
    "\n"
    "timeFormat         general;\n"
    "\n"
    "writeCompression   off;\n"
    "\n"
    "timeFormat         general;\n"
    "\n"
    "timePrecision      6;\n"
    "\n"
    "runTimeModifiable  true;\n"
    "\n"
    "\n";


}
void blockMeshGen::generateSurfaceFeature()
{
    out.close();
    out.open("boundary/surfaceFeatureExtractDict");

    out << "\n"
    "FoamFile\n"
    "{\n"
    "   version     2.0;\n"
    "   format      ascii;\n"
    "   class       dictionary;\n"
    "   object      surfaceFeatureExtractDict;\n"
    "}\n\n"
    "blade.stl\n"  
    "{\n"  
    "   extractionMethod    extractFromSurface;\n"  
    "   extractFromSurfaceCoeffs\n"  
    "   {\n"  
    "       includedAngle 180;\n"  
    "   }\n"  
    "\n"  
    "   subsetFeatures\n"  
    "   {\n"  
    "       nonManifoldEdges    yes;\n"  
    "       openEdges           yes;\n"  
    "   }\n"  
    "   writeObj    yes;\n"  
    "}\n"; 
}
void blockMeshGen::generateCreatePatch()
{
    std::ofstream output("systemFile/createPatchDict");

    output << "FoamFile\n"
    "{\n"
    "   version     2.0;\n"
    "   format      ascii;\n"
    "   class       dictionary;\n"
    "   object      createPatchDict;\n"
    "}\n\n"
    "\n"
    "pointSync      false;\n" 
    "\n" 
    "patches\n" 
    "(\n" 
    "   {\n" 
    "       name            cyclicFace1;\n" 
    "       constructFrom   patches;\n" 
    "\n" 
    "       patchInfo\n" 
    "       {\n" 
    "           type            cyclic;\n" 
    "           matchTolerance  0.01;\n" 
    "           neighbourPatch  cyclicFace2;\n" 
    "           transform       rotational;\n" 
    "           rotationAxis    (0 0 1);\n" 
    "           rotationCentre  (0 0 0);\n" 
    "\n" 
    "       }\n" 
    "\n" 
    "   patches (cyclic1);\n" 
    "\n" 
    "   set         cyclicFace1;\n" 
    "   }\n" 
    "\n" 
    "\n" 
    "\n" 
    "\n" 
    "\n" 
        "   {\n" 
    "       name            cyclicFace2;\n" 
    "       constructFrom   patches;\n" 
    "\n" 
    "       patchInfo\n" 
    "       {\n" 
    "           type            cyclic;\n" 
    "           matchTolerance  0.01;\n" 
    "           neighbourPatch  cyclicFace1;\n" 
    "           transform       rotational;\n" 
    "           rotationAxis    (0 0 1);\n" 
    "           rotationCentre  (0 0 0);\n" 
    "\n" 
    "       }\n" 
    "\n" 
    "   patches (cyclic2);\n" 
    "\n" 
    "   set         cyclicFace2;\n" 
    "   }\n" "\n" 
    "\n" 
    ");\n" 
    "\n" 
    "\n" 
    "\n" 
    "\n" 
    "\n" 
    "\n"; 


}

void blockMeshGen::generateObj(int resolution)
{
std::ofstream obj("stl/blade.obj");

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

                    obj << "f " << d << " " << d+1 << " " << d+resolution+1 << "\n";

                }

                if(l == 1)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[d+resolution+ 2]);
                    tempValue.push_back(vertices[d+1]); 
                    tempValue.push_back(vertices[d+resolution+1]);

                    obj << "f " << d+resolution+1 << " " << d+1 << " " << d+resolution+2 << "\n";

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
            
            printCoordinatesObj(obj);
         
        
        }
    
    }

    int limitB = vertices.size() - resolution -1;
    for(int d = 0; d < resolution - 0; d++)
    {
        for(int l = 0; l < 2; l++)
        {

                if(l == 0 && remainder(d + 0 , resolution/2) != 0)
                {   
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB + d]);
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d+1]);

                    obj << "f " << limitB - d + 1<< " " << limitB + d << " " << limitB + d+resolution << "\n";

                }

                if(l == 1 && remainder(d + 1 , resolution/2) != 0)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB + d+1]); 
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB - d + resolution - 1]);

                    obj << "f " << limitB - d + resolution - 1 << " " << limitB - d + resolution << " " << limitB + d+1 << "\n";

                }
            

            printCoordinatesObj(obj);
         
        }
    
    }

    limitB = 0;
    for(int d = 1; d < resolution - 1; d++)
    {

        for(int l = 0; l < 2; l++)
        {

                if(l == 0 && remainder(d + 0 , resolution/2) != 0)
                {   
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d]);
                    tempValue.push_back(vertices[limitB + d+1]);

                    obj << "f " << limitB - d + 1 << " " << limitB + d << " " << limitB + d+resolution << "\n";

                }

                if(l == 1 && remainder(d + 1 , resolution/2) != 0)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d+1]); 
                    tempValue.push_back(vertices[limitB - d + resolution - 1]);

                    obj << "f " << limitB - d + resolution - 1 << " " << limitB - d + 1 << " " << limitB + d+resolution << "\n";
                    
                }
            

            printCoordinatesObj(obj);
         
        }
    
    }

}

void blockMeshGen::printCoordinatesObj(std::ofstream& obj)
{
    for(int i = 0; i < 3; i++)
    {
    obj << "v " << tempValue[i][0] << " " << tempValue[i][1] << " " << tempValue[i][2] << "\n";
    }
}

void blockMeshGen::getSeparationVector(double input)
{
    separationVector = input;
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
            
            printCoordinates();
         
        
        }
    
    }

    int limitB = vertices.size() - resolution -1;
    for(int d = 0; d < resolution - 0; d++)
    {
        for(int l = 0; l < 2; l++)
        {

                if(l == 0 && remainder(d + 0 , resolution/2) != 0)
                {   
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB + d]);
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d+1]);
                }

                if(l == 1 && remainder(d + 1 , resolution/2) != 0)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB + d+1]); 
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB - d + resolution - 1]);
                }
            

            printCoordinates();
         
        }
    
    }

    limitB = 0;
    for(int d = 1; d < resolution - 1; d++)
    {

        for(int l = 0; l < 2; l++)
        {

                if(l == 0 && remainder(d + 0 , resolution/2) != 0)
                {   
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d]);
                    tempValue.push_back(vertices[limitB + d+1]);
                }

                if(l == 1 && remainder(d + 1 , resolution/2) != 0)
                {
                    tempValue.clear();
                    tempValue.push_back(vertices[limitB - d + resolution]);
                    tempValue.push_back(vertices[limitB + d+1]); 
                    tempValue.push_back(vertices[limitB - d + resolution - 1]);
                }
            

            printCoordinates();
         
        }
    
    }

    out << "endsolid blade\n";

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
    out.open("stl/inlet.stl");

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
    out.open("stl/outlet.stl");

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
    out.open("stl/bot.stl");

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
    out.open("stl/top.stl");

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