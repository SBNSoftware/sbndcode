////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesPCA.cpp
//
// \brief Definition of TPCLinesPCA
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPCLINES_PCA_H
#define TPCLINES_PCA_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"

class TPCLinesPCA {
    public:
        TPCLinesPCA() {}

        LineEquation PerformPCA2D(std::vector<SHit>& data);

        LineEquation PerformPCA2DThreshold(std::vector<SHit> &hitList, float th, bool reverse = false);
};

#endif // TPCLINES_PCA_H