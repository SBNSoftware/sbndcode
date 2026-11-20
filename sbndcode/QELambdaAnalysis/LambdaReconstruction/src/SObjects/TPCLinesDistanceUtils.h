////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef DISTANCE_UTILS_H
#define DISTANCE_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>


#include "TPCSimpleClusters.h"
#include "TPCLinesHitDistanceUtils.h"

namespace TPCLinesDistanceUtils{

    double GetClusterMinDistance(SCluster cluster1, SCluster cluster2);

    double GetClusterMinDistanceX(SCluster cluster1, SCluster cluster2);

    double GetClusterConnectedness(SCluster cluster1, SCluster cluster2);

    double GetClusterConnectednessOverlap(SCluster cluster1, SCluster cluster2);

    double GetpClusterDistanceW(SHit p1, SHit p2);

    double GetpClusterDistanceOverlap(SHit p1, SHit p2);
}

#endif // TPC_SIMPLE_CLUSTERS_H