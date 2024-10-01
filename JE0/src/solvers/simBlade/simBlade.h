#ifndef SIMBLADE_H
#define SIMBLADE_H

#include "../aeroBlade/aeroBlade.h"

namespace simBlade
{
    void storeInDatabaseRecursive();
    void generateAerofoilModel(int i, int j);
};

#endif