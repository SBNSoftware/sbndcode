////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.h
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_CLUSTERS_H
#define TPC_SIMPLE_CLUSTERS_H

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iosfwd>
#include <limits>

#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"
#include "TPCLinesPCA.h"
#include "TPCLinesHitDistanceUtils.h"


class SCluster {
    private:
        int fNHits;

        std::vector<SHit> fHitList;
        std::vector<double> fConnectednessV;
        std::vector<double> fCompactnessV;
        std::vector<double> fConnectedness1DV;
        std::vector<double> fWidths;
        
        int fMainClusterId;
        std::vector<int> fClusterIDs;

        double fCompactness;
        double fCompactnessRMS;
        double fConnectedness;
        double fConnectednessRMS;
        double fConnectedness1D;
        double fConnectedness1DRMS;
        double fAverageWidth;

        static constexpr double DefaultMaxDistance = std::numeric_limits<double>::max();
    
    public:
        SCluster(std::vector<SHit> hitList={});

        int NHits(){return fNHits;}
        void AddHit(SHit& hit);
        
        // return methods
        std::vector<SHit> GetHits(){return fHitList;};
        std::vector<double> GetConnectednessV(){return fConnectednessV;};
        double GetConnectedness(){return fConnectedness;};
        double GetConnectednessRMS(){return fConnectednessRMS;};
        std::vector<double> GetCompactnessV(){return fCompactnessV;};
        double GetCompactness(){return fCompactness;};
        double GetCompactnessRMS(){return fCompactnessRMS;};
        double GetAverageWidth(){return fAverageWidth;};
        int GetMainClusterID() {return fMainClusterId;};
        void CalculateMainClusterID();
        std::vector<int> GetClusterIDs() {return fClusterIDs;};

        // get hit distance to clusters method
        template <typename T> double GetMinDistanceToCluster(const T& h);
        template <typename T> double GetMinDistanceToClusterW(const T& h);
        template <typename T> double GetMinDistanceToClusterOverlap(const T& h);
        template <typename T> double GetMinDistanceToCluster1D(const T& h);
        template <typename T> std::pair<SHit, double>  GetClosestHitToPoint(const T& h);
        

        friend std::ostream& operator<<(std::ostream& out, SCluster & c);

};




// Class to represent a linear cluster
class SLinearCluster {
    private:
        int fId;
        SCluster fHitCluster;

        static constexpr double DefaultMax = std::numeric_limits<int>::max();
   
        float fMinX;
        float fMinY;
        float fMaxX;
        float fMaxY;
        float fYAtMaxX;
        float fYAtMinX;
        float fMeanX;
        float fMeanY;

        LineEquation fTrackEquation;
        LineEquation fTrackEquationStart;
        LineEquation fTrackEquationEnd;
        std::map<int, LineEquation> fTrackEquationSlidingWindowV;
        
        SCluster fResidualHitCluster;
        SCluster fMainHitCluster;

        SPoint fStartPoint;
        SPoint fEndPoint;
        SPoint fCoMPoint;
        bool fHasStartEndPoints;
        bool fHasResidualHits;

        double CalculateSlope(const SHit& p1, const SHit& p2);


    public:
        
        // constructor
        SLinearCluster(std::vector<SHit> hitList={});

        // return methods
        int GetId() {return fId;}
        std::vector<SHit> GetHits() {return fHitCluster.GetHits();}
        int NHits() {return fHitCluster.NHits();}
        double GetIntegral();

        // Getter functions
        float GetMinX()  { return fMinX; }
        float GetMinY()  { return fMinY; }
        float GetMaxX()  { return fMaxX; }
        float GetMaxY()  { return fMaxY; }
        float GetMeanX()  { return fMeanX; }
        float GetMeanY()  { return fMeanY; }
        float GetYatMinX()  { return fYAtMinX; }
        float GetYatMaxX()  { return fYAtMaxX; }
        SPoint GetStartPoint()  { return fStartPoint; }
        SPoint GetEndPoint()  { return fEndPoint; }
        SPoint GetCoMPoint()  { return fCoMPoint; }

        float GetHitDensity();

        LineEquation GetTrackEquation() { return fTrackEquation; }
        LineEquation GetTrackEquationStart() { return fTrackEquationStart; }
        LineEquation GetTrackEquationEnd() { return fTrackEquationEnd; }

        SHit FindMaxVariationPointSlidingWindow(std::vector<SHit> data, size_t N);

        SCluster GetHitCluster(){return fHitCluster;};
    
        std::vector<double> GetConnectednessV(){return fHitCluster.GetConnectednessV();};
        double GetConnectedness(){return fHitCluster.GetConnectedness();};
        std::vector<double> GetCompactnessV(){return fHitCluster.GetCompactnessV();};
        double GetCompactness(){return fHitCluster.GetCompactness();};
        double GetAverageWidth(){return fHitCluster.GetAverageWidth();};
    
        void AssignId(int id) {fId = id;}
        void AddHit(SHit hit) {fHitCluster.AddHit(hit);}

        std::vector<int> detect_outliers_iqr(std::vector<float> data, float threshold);
        std::vector<int> detect_outliers_iqr2(std::vector<double> data, double threshold = 1.5);

        void FillResidualHits(bool customKinkPoint=false);

        double GetOccupancy();
        double GetOccupancy1D();

        SPoint GetEdgeHit(SPoint p);
        SHit GetOppositeEdgeHit(SPoint p);

        double GetLength();

        void FillSlidingWindowLineEquations(size_t N);
        LineEquation GetLineEquationAtX(int x);
};


#endif // TPC_SIMPLE_CLUSTERS_H