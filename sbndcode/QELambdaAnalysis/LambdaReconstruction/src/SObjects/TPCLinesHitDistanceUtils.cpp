////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesHitDistanceUtils.h"

namespace TPCLinesDistanceUtils{
    // Mean and Std Dev function for std::vector
    double CalculateMean( std::vector<double>& values) {

        double sum = 0;
        for ( double& value : values) {
            sum += value;
        }
        return sum / values.size();
    }

    double CalculateStdDev( std::vector<double>& values) {

        double mean = CalculateMean(values);
        double sumSqDiff = 0;
        for ( double& value : values) {
            double diff = value - mean;
            sumSqDiff += diff * diff;
        }
        return std::sqrt(sumSqDiff / values.size());
    }

    template <typename T1, typename T2>
    double GetHitDistance(const T1& hit1, const T2& hit2) {
        double xdif = hit1.X() - hit2.X();
        double ydif = hit1.Y() - hit2.Y();

        return std::sqrt(xdif * xdif + ydif * ydif);
    }


    // Hit distance funcitons
    double GetPointDistance( SPoint p1,  SPoint p2) {

        double xdif = p1.X() - p2.X();
        double ydif = p1.Y() - p2.Y();

        return std::sqrt( xdif*xdif + ydif*ydif);
    }

    template <typename T1, typename T2>
    int GetHitDistance1D(const T1& hit1, const T2& hit2) {
        return std::abs(hit1.X() - hit2.X());
    }


    template <typename T1, typename T2>
    double GetHitDistanceW(const T1& hit1, const T2& hit2) {
        
        // distance in X
        double dX = std::pow(hit1.X() - hit2.X(), 2);
        
        // distance between centers
        double d0 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y(), 2));
        
        // dist of hit 1 to width of hit 2
        double dp1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit2.Width(), 2));
        double dm1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit2.Width(), 2));
        double d1 = std::min(dp1, dm1);
        
        // dist of hit 2  to width of hit 1
        double dp2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit1.Width(), 2));
        double dm2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit1.Width(), 2));
        double d2 = std::min(dp2, dm2);
        
        // return min
        return std::min(d0, (d1 + d2) / 2.);
    }
    

    template <typename T1, typename T2>
    double GetHitDistanceOverlap(const T1& hit1, const T2& hit2) {
        double y_range1_min = hit1.Y() - hit1.Width();
        double y_range1_max = hit1.Y() + hit1.Width();
        double y_range2_min = hit2.Y() - hit2.Width();
        double y_range2_max = hit2.Y() + hit2.Width();
        bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
        double dX = ( std::abs(hit1.X() - hit2.X())<=1 )? 0:std::pow(hit1.X() - hit2.X(), 2);
        double dY = 0;
        if (!overlap) {
            double d1 = std::abs(y_range1_min-y_range2_max);
            double d2 = std::abs(y_range2_min-y_range1_max);
            dY = std::min(d1, d2);
        }
        return std::sqrt(dX + dY);
    }

    


    bool HitWidthOverlap(SHit& hit1, SHit& hit2) {
        std::pair<double, double> y_range1 = std::make_pair( hit1.Y() - hit1.Width(), hit1.Y() + hit1.Width() );
        std::pair<double, double> y_range2 = std::make_pair( hit2.Y() - hit2.Width(), hit2.Y() + hit2.Width() );
        
        return y_range1.first <= y_range2.second && y_range1.second >= y_range2.first;
    }


    // instances
    template double GetHitDistance(const SHit& hit1, const SHit& hit2);
    template double GetHitDistance(const SPoint& hit1, const SPoint& hit2);
    template double GetHitDistance(const SHit& hit1, const SPoint& hit2);

    template int GetHitDistance1D(const SHit& hit1, const SHit& hit2);
    template int GetHitDistance1D(const SPoint& hit1, const SPoint& hit2);
    template int GetHitDistance1D(const SHit& hit1, const SPoint& hit2);
    
    template double GetHitDistanceW(const SHit& hit1, const SHit& hit2);
    template double GetHitDistanceOverlap(const SHit& hit1, const SHit& hit2);

}