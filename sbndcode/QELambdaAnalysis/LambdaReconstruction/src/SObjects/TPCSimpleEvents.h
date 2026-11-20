////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleEvents.h
//
// \brief Definition of SimpleEvent
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_EVENT_H
#define TPC_SIMPLE_EVENT_H

#include <iostream>
#include <string>
#include <iosfwd>
#include <unordered_set>



#include "TPCSimpleLines.h"
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"



class SOrigin {
    private:
        SPoint fVertex;
        double fYError;
        std::vector<SLinearCluster> fTrackList;
        int fMultiplicity;
        int fNHits;
        bool fEdgeOrigin;
        std::unordered_set<int> fTrackConnectionIDs;
              
    public:
        SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge, double yError, int kinkParentId = -1);
        
        SPoint GetPoint(){ return fVertex;};
        double GetYError(){ return fYError;};
        std::vector<SLinearCluster> GetTracks(){ return fTrackList; };
        SLinearCluster GetTrack(int ix){ return fTrackList[ix];};
        std::vector<int> GetTracksIDs();
        std::vector<int> GetTracksConnectionsIDs();
        int Multiplicity(){return fMultiplicity;};
        void AddTrack(SLinearCluster track, SPoint p, double yError); 
        bool IsEdgeOrigin(){return fEdgeOrigin;};
        bool HasTrackIndex(int ix);
        bool IsConnectedToTrackIndex(int ix);
        int NHits(){return fNHits;};
        double TotalCharge();
        SLinearCluster GetTrackEntry(size_t ix){ return fTrackList[ix];};

        //friend std::ostream& operator<<(std::ostream& out, SOrigin const& ori);

        friend std::ostream& operator<<(std::ostream& out, SOrigin & ori)
        {
            out << " SOrigin -- Vertex=(" << ori.GetPoint().X() << ", " <<ori.GetPoint().Y() << ")"<< " Mult="<<ori.Multiplicity()<< " Tracks= ";
            for(auto & trk: ori.GetTracks()){
                out<<trk.GetId()<<" ";
            }
            out<<"  NHits"<<ori.NHits()<<std::endl;
            return out;
        }
};



class SEvent {
    private:
        std::vector<SLinearCluster> fTrackList;
        std::vector<SOrigin> fOriginList;
        std::vector<STriangle> fAngleList;
        std::vector<SOrigin> fAssociatedOrigins;
        double fHitDensity;
        std::vector<SHit> fFreeHits;

        std::map<int, std::vector<int>> fTrackConnections;

    public:
        SEvent(std::vector<SLinearCluster> tracks={}, std::vector<SOrigin> origins={}, std::vector<STriangle> angles = {}, std::vector<SOrigin> associatedOrigins = {}, double hitDensity=0, std::vector<SHit> freeHits = {});

        std::vector<SOrigin> GetOrigins(){return fOriginList;};
        int GetNOrigins();
        int GetNOriginsMult(int mult);
        int GetNOriginsMultGt(int mult);
        int GetNOriginsMultGt(int mult, int id1, int id2);

        int GetNAngles(){ return fAngleList.size();};

        std::vector<STriangle> GetAngleList(){return fAngleList;};
        std::vector<SOrigin> GetAssociatedOrigins(){return fAssociatedOrigins;};

        double HitDensity(){return fHitDensity;};

        void PrintTrackConnections();

        int NHits();
        int NFreeHits(){return fFreeHits.size();};

        bool IsOriginAssociatedToTrack(SOrigin ori, int trackId);
        bool IsTrackAssociatedToTrack(int trackId, int trackId2);

        // Function to return the free hits
        std::vector<SHit> GetFreeHits() {return fFreeHits;};
        
        // Function to return the tracks
        std::vector<SLinearCluster> GetTracks(){return fTrackList;};
        
        SEvent UnassociatedOrigins(STriangle triangle);

        void FreeHitsAroundTriangle(STriangle triangle,
                                    int & nDirtHitsInTriangle,
                                    double & nFractionDirtHitsInTriangle,
                                    int & nDirtHitsInTriangleWires,
                                    double & nFractionDirtHitsInTriangleWires);

        void GetUnassociatedHits(STriangle triangle, int &nFreeHits, int &fNUnassociatedHits);                            
};

#endif // TPC_SIMPLE_HITS_H