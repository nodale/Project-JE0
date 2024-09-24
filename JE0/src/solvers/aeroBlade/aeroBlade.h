#ifndef AEROBLADE_H
#define AEROBLADE_H

#include "../infoBlade/infoBlade.h"
#include "../thermoBlade/thermoBlade.h"

namespace aeroBlade
{

extern dVec<double> I, J, K, L;
extern sVec<double> Vn, Vt, phi, lambda, Sj;

extern sVec<double> pointX, pointY, midPointX, midPointY;

extern double aeroGamma; 

void init();

};

#endif