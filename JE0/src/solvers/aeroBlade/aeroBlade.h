#ifndef AEROBLADE_H
#define AEROBLADE_H

#include "../infoBlade/infoBlade.h"
#include "../thermoBlade/thermoBlade.h"

namespace aeroBlade
{

extern std::vector<std::vector<double>> I, J, K, L;
extern std::vector<double> Vn, Vt, phi, lambda, Sj;

extern std::vector<double> pointX, pointY;

void init();

};

#endif