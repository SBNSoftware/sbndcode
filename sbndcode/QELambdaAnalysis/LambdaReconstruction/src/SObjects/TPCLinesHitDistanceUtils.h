////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef HITDISTANCE_UTILS_H
#define HITDISTANCE_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "TPCSimpleHits.h"

namespace TPCLinesDistanceUtils{
    // Mean and Std Dev function for std::vector
    double CalculateMean( std::vector<double>& values);

    double CalculateStdDev( std::vector<double>& values);

    template <typename T1, typename T2>
    double GetHitDistance(const T1& hit1, const T2& hit2);

    // Hit distance funcitons
    double GetPointDistance( SPoint p1,  SPoint p2);

    template <typename T1, typename T2>
    int GetHitDistance1D(const T1& hit1, const T2& hit2);

    template <typename T1, typename T2>
    double GetHitDistanceW(const T1& hit1, const T2& hit2);

    template <typename T1, typename T2>
    double GetHitDistanceOverlap(const T1& hit1, const T2& hit2);

    bool HitWidthOverlap(SHit& hit1, SHit& hit2);
}

#endif // TPC_SIMPLE_CLUSTERS_H