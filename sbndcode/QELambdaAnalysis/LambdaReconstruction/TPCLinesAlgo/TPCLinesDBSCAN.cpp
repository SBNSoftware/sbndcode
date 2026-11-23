////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesDBSCAN.cpp
//
// \brief Definition of TPCLinesDBSCAN
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesDBSCAN.h"

DBSCAN::DBSCAN(double epsilon, int minPts) 
    :epsilon(epsilon),
    minPts(minPts)
{}

void DBSCAN::setDistanceFunction(double (*distFunc)( SHit&,  SHit&)) {
    distanceFunction = distFunc;
}


std::vector<int>& DBSCAN::getClusterAssignment()  {
    return clusterAssignment;
}

void DBSCAN::fit( std::vector<SHit>& points) {
    clusterAssignment.assign(points.size(), 0);  // 0 indicates unassigned
    
    int clusterIdx = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (clusterAssignment[i] == 0) {  // Not assigned to any cluster
            std::vector<int> neighbors = regionQuery(points, points[i]);
            if ((int)neighbors.size() < minPts) {
                clusterAssignment[i] = -1;  // Mark as noise
            }
            else {
                ++clusterIdx;
                expandCluster(points, i, clusterIdx);
            }
        }
    }
}

std::vector<int> DBSCAN::regionQuery( std::vector<SHit>& points,  SHit& p) {
    std::vector<int> neighbors;
    for (size_t i = 0; i < points.size(); ++i) {
        if (distanceFunction(p, points[i]) <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void DBSCAN::expandCluster( std::vector<SHit>& points, int pointIdx, int clusterIdx) {
   
    
    std::vector<int> seeds = regionQuery(points, points[pointIdx]);

    clusterAssignment[pointIdx] = clusterIdx;
    
    for (size_t i = 0; i < seeds.size(); ++i) {
        int currentIdx = seeds[i];
        // not assigned to any cluster
        if (clusterAssignment[currentIdx] == 0) {
            clusterAssignment[currentIdx] = clusterIdx;
            std::vector<int> currentSHitNeighbors = regionQuery(points, points[currentIdx]);
            if ((int)currentSHitNeighbors.size() >= minPts) {
                expandCluster(points, currentIdx, clusterIdx);
            }
        }
        // noise point
        if (clusterAssignment[currentIdx] == -1) {
            clusterAssignment[currentIdx] = clusterIdx;
        }
    }
}

// --------------- Distance functions definitions for DBSCAN

double DBSCANHitEuclidianDistance( SHit& p1,  SHit& p2) {
    //std::cout<<" In euclidian distance\n";
    return std::sqrt((p1.X() - p2.X()) * (p1.X() - p2.X()) + (p1.Y() - p2.Y()) * (p1.Y() - p2.Y()));
}

double DBSCANHitManhattanDistance( SHit& p1,  SHit& p2) {
    return std::abs(p1.X() - p2.X()) + std::abs(p1.Y() - p2.Y());
}

double DBSCANHitWidthDistance( SHit& p1, SHit& p2) {
    //std::cout<<" In width distance\n";
    double dX = std::pow(p1.X() - p2.X(), 2);
    
    double d0 = std::sqrt( dX + std::pow(p1.Y() - p2.Y(), 2));
    
    double dY1 =  std::min( std::abs(p1.Y() + p1.Width() - p2.Y()),  std::abs(p1.Y() - p1.Width() - p2.Y()) );
    double dY2 =  std::min( std::abs(p1.Y() - p2.Y() + p2.Width()),  std::abs(p1.Y() - p2.Y() - p2.Width()) );

    double d1 = std::sqrt( dX + std::pow(dY1, 2));
    double d2 = std::sqrt( dX + std::pow(dY2, 2));
    
    //return std::min(d0, std::min( d1, d2 ));
    return std::min(d0, 0.5*(d1+d2));
}

double DBSCANHitWidthDistance2( SHit& p1, SHit& p2) {
    //std::cout<<" In width distance\n";
    double d0 = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
    double dp = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() + p2.Width(), 2));
    double dm = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y() - p2.Width(), 2));
    return std::min(d0, std::min(dp, dm));
}

double DBSCANHitEuclidianDistanceDriftConversion( SHit& p1,  SHit& p2) {
    //std::cout<<" In euclidian distance\n";
    return std::sqrt((p1.X() - p2.X()) * (p1.X() - p2.X()) + 0.075*(p1.Y() - p2.Y()) * (p1.Y() - p2.Y()));
}



double DBSCANHitOverlapDistance( SHit& p1, SHit& p2) {
    double y_range1_min = p1.Y() - p1.Width();
    double y_range1_max = p1.Y() + p1.Width();
    double y_range2_min = p2.Y() - p2.Width();
    double y_range2_max = p2.Y() + p2.Width();
    bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
    double dX = p1.X() - p2.X();
    double dY = 0;

    if (!overlap) {
        dY = p1.Y() - p2.Y();
    }
    
    return std::sqrt(std::pow(dX, 2) + std::pow(dY, 2));
}
