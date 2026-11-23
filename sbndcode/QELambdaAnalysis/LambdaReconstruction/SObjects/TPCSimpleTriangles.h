////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleTriangle.h
//
// \brief Definition of SimpleTriangle
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_TRIANGLE_H
#define TPC_SIMPLE_TRIANGLE_H

#include <iostream>
#include <string>
#include <iosfwd>
#include <map>

#include "TPCSimpleLines.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleHits.h"
#include "TPCLinesHitDistanceUtils.h"


class STriangle {
    private:
        SPoint fMainVertex;
        SHit fMainVertexHit;
        SPoint fVertexB;
        SPoint fVertexC;
        SPoint fMidPoint;
        SPoint fDirectorVector;
        LineEquation fDirection;
        LineEquation fDir1;
        LineEquation fDir2;
        LineEquation fMomentumHypo1;
        LineEquation fMomentumHypo2;

        SLinearCluster fTrack1;
        SLinearCluster fTrack2;
        SLinearCluster fMainTrack;

        double fOpeningAngle;
        
    public:
        STriangle();

        STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, SLinearCluster track1, SLinearCluster track2, SLinearCluster mainTrack, double weight_b=1, double weight_c=1);

        SLinearCluster GetTrack1() const {
            return fTrack1;
        }

        SLinearCluster GetTrack2() const {
            return fTrack2;
        }

        SLinearCluster GetMainTrack() const {
            return fMainTrack;
        }
        
        SPoint GetMainVertex() const {
            return fMainVertex;
        }

        SHit GetMainVertexHit() const {
            return fMainVertexHit;
        }

        SPoint GetVertexB() const {
            return fVertexB;
        }

        SPoint GetVertexC() const {
            return fVertexC;
        }

        SPoint GetMidPoint() const {
            return fMidPoint;
        }

        SPoint GetDirectorVector() const {
            return fDirectorVector;
        }

        LineEquation GetDirection() const {
            return fDirection;
        }

        LineEquation GetMomentumHypo1() const {
            return fMomentumHypo1;
        }

        LineEquation GetMomentumHypo2() const {
            return fMomentumHypo2;
        }

        double GetOpeningAngle() const {
            return fOpeningAngle;
        }

        double GetSideLenghtB() const {
            return std::hypot( fMainVertex.X() - fVertexB.X(), fMainVertex.Y() - fVertexB.Y() );
        }

        double GetSideLenghtC() const {
            return std::hypot( fMainVertex.X() - fVertexC.X(), fMainVertex.Y() - fVertexC.Y() );
        }

        double GetOppositeSideLenght() const {
            return std::hypot( fVertexB.X() - fVertexC.X(), fVertexB.Y() - fVertexC.Y() );
        }
        
        // Get NHits functions
        int GetNHitsTriangle();
        int GetNHitsTrack1();
        int GetNHitsTrack2();
        int GetNHitsMainTrack();

        // Get all hits function
        std::vector<SHit> GetAllHits();
        std::vector<SHit> GetMainTrackAllHits();

        // Get Length functions
        double GetLengthTrack1();
        double GetLengthTrack2();
        double GetLengthMainTrack();

        double GetDecayAngleDifference();

        int GetNWires();

        double ComputeCoveredArea();

        double GetTriangleMAE();

        double GetGap();

        double GetMinX();
        double GetMaxX();
        double GetMinY();
        double GetMaxY();

        double GetOverlapWithMainTrack();
        
        std::vector<int> GetClusterIDs();

        std::vector<int> GetMainTrackClusterIDs();
        std::vector<int> GetTrack1ClusterIDs();
        std::vector<int> GetTrack2ClusterIDs();


        void GetVertexXYZ(double &x, double &y, double &z);
        void GetPrimaryLeptonVertexXYZ(double &x, double &y, double &z);
        
};

#endif // TPC_SIMPLE_HITS_H