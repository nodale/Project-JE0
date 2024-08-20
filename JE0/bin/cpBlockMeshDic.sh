#cp boundary/blockMeshDict ~/Documents/test_ground/system/
mkdir ~/Documents/JE0/CFD_cases/snappyHexMesh_test/constant
mkdir ~/Documents/JE0/CFD_cases/snappyHexMesh_test/constant/polyMesh/
mkdir ~/Documents/JE0/CFD_cases/snappyHexMesh_test/constant/triSurface/
mkdir ~/Documents/JE0/CFD_cases/snappyHexMesh_test/system
#cp boundary/boundary ~/Documents/test_ground/constant/polyMesh/
cp -a stl/. ~/Documents/JE0/CFD_cases/snappyHexMesh_test/constant/triSurface/
#cp CPD/snappyHexMeshDict ~/Documents/test_ground/system/
#cp boundary/surfaceFeatureExtractDict ~/Documents/test_ground/system/
#cp boundary/createPatchDict ~/Documents/test_ground/system/
#cp boundary/controlDict ~/Documents/test_ground/system/
cp -a systemFile/. ~/Documents/JE0/CFD_cases/snappyHexMesh_test/system/
