#include "Cl_Cd.h"
	
template<typename d>
CL_CD::CL_CD(d a, d b)
: alpha {a}
, beta {b}
{
	gamma = beta - alpha;
	//Cl's range is between 0.1 to 0.8
	Cl = (0.1077 * gamma) + 0.095;
	//Cd's range is between 0.0 to 0.07
	Cd = 0.007 * gamma * gamma + 0.005 * gamma + 0.005;
};

template<typename d>
d CL_CD::getCl(d a, d b)
{
	gamma = beta - alpha;

	Cl = (0.1077 * gamma) + 0.095;

	return Cl;
};

template<typename d>
d CL_CD::getCd(d a, d b)
{
	gamma = beta - alpha;

	Cd = 0.007 * gamma * gamma + 0.005 * gamma + 0.005;

	return Cd;
};
