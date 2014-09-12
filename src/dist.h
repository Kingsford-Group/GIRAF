#ifndef DIST_H
#define DIST_H

#include <fstream>
#include "tree.h"

void ComputePairDistances(istream &, istream &, bool, ostream *, DistanceMatrix &,
        DistanceMatrix &);

#endif
