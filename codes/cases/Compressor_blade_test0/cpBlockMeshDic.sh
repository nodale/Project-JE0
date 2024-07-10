#cp boundary/blockMeshDict ~/Documents/test_ground/system/
mkdir ~/Documents/test_ground/constant
mkdir ~/Documents/test_ground/constant/polyMesh/
mkdir ~/Documents/test_ground/constant/triSurface/
mkdir ~/Documents/test_ground/system
#cp boundary/boundary ~/Documents/test_ground/constant/polyMesh/
cp -a stl/. ~/Documents/test_ground/constant/triSurface/
#cp CPD/snappyHexMeshDict ~/Documents/test_ground/system/
#cp boundary/surfaceFeatureExtractDict ~/Documents/test_ground/system/
#cp boundary/createPatchDict ~/Documents/test_ground/system/
#cp boundary/controlDict ~/Documents/test_ground/system/
cp -a systemFile/. ~/Documents/test_ground/system/
