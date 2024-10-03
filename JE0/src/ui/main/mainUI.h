#ifndef MAINUI_H
#define MAINUI_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../solvers/infoBlade/infoBlade.h"
#include "../../solvers/thermoBlade/thermoBlade.h"
#include "../../solvers/simBlade/simBlade.h"

void quit();
void setupDatabase();
void initialiseInitCond();
void initialiseThermoVar();
void initialiseThermoVar_v2();
void configureCompressionRatio(std::istringstream& stream);
void configureInitAlpha(std::istringstream& stream);
void deleteDatabase();

#endif