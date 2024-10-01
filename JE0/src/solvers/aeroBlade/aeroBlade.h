#ifndef AEROBLADE_H
#define AEROBLADE_H

#include "../infoBlade/infoBlade.h"
#include "../thermoBlade/thermoBlade.h"

namespace aeroBlade
{

extern dVec<double> I, J, K, L;
extern sVec<double> Vn, Vt, phi, lambda, Sj;

extern sVec<double> pointX, pointY, midPointX, midPointY;

extern double aeroGamma, aeroCl;

void sourceVortexPanelMethod(int i, int j, int r);
void genBlade(int i, int j);
void findCombinationAlpha(int sampleSize, int maxTries);
void drawDeHallersNumber();
void drawDeHallersNumber_v2();
std::complex<double> joukowskyTransform(std::complex<double> z, double shape, std::complex<double> thetaC);
void getAoA(int i, int j, int R);
void storeInDatabaseRecursive();
void drawEverthing(int I, int j, int R);

};

#endif