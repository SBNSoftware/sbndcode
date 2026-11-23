////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesDistanceUtils.h"

namespace TPCLinesDistanceUtils{

    double GetClusterMinDistance(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = std::sqrt(std::pow(hit1.X() - hit2.X(), 2) + std::pow(hit1.Y() - hit2.Y(), 2));
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }

    double GetClusterMinDistanceX(SCluster cluster1, SCluster cluster2) {
        int minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = std::abs(hit1.X() - hit2.X());
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }

    double GetClusterConnectedness(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = GetHitDistanceW(hit1, hit2);
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }

    double GetClusterConnectednessOverlap(SCluster cluster1, SCluster cluster2) {
        double minDistance = 1e4;
        for (SHit& hit1 : cluster1.GetHits()) {
            for (SHit& hit2 : cluster2.GetHits()) {
                double distance = GetHitDistanceOverlap(hit1, hit2);
                if (distance < minDistance) {
                    minDistance = distance;
                }
            }
        }
        return minDistance;
    }
    


    double GetpClusterDistanceW(SHit p1, SHit p2) {
        double dX = std::pow(p1.X() - p2.X(), 2);
    
        double d0 = std::sqrt( dX + std::pow(p1.Y() - p2.Y(), 2));
        
        double dY1 =  std::min( std::abs(p1.Y() + p1.Width() - p2.Y()),  std::abs(p1.Y() - p1.Width() - p2.Y()) );
        double dY2 =  std::min( std::abs(p1.Y() - p2.Y() + p2.Width()),  std::abs(p1.Y() - p2.Y() - p2.Width()) );

        double d1 = std::sqrt( dX + std::pow(dY1, 2));
        double d2 = std::sqrt( dX + std::pow(dY2, 2));
        
        return std::min(d0, std::min( d1, d2 ));
        //return std::min(d0, 0.5*(d1+d2));
        //return std::min(d0, dY2);
    }


    double GetpClusterDistanceOverlap(SHit p1, SHit p2) {
        double y_range1_min = p1.Y() - p1.Width();
        double y_range1_max = p1.Y() + p1.Width();
        double y_range2_min = p2.Y() - p2.Width();
        double y_range2_max = p2.Y() + p2.Width();
        bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
        double dX = p1.X() - p2.X();
        double dY = 1;
        if (!overlap) {
            dY = std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
        }
        return std::sqrt(std::pow(dX, 2) + std::pow(dY, 2));
    }

}
