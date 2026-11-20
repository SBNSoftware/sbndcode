////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleEvents.cpp
//
// \brief Definition of SimpleEvents
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleEvents.h"

SOrigin::SOrigin(SPoint p, std::vector<SLinearCluster> tracks, bool isEdge, double yError, int kinkParentId) :
    fVertex(p),
    fYError(yError),
    fTrackList(tracks),
    fMultiplicity(tracks.size()),
    fNHits(0),
    fEdgeOrigin(isEdge)
{

    for(SLinearCluster & trk:tracks){
        fTrackConnectionIDs.insert(trk.GetId());
    }
    if(kinkParentId!=-1){
        fTrackConnectionIDs.insert(kinkParentId);
    }
    
    for(SLinearCluster & trk:tracks){
        fNHits+=trk.NHits();
    }
}

std::vector<int> SOrigin::GetTracksIDs(){
    std::vector<int> ids;
    for(const int &id: fTrackConnectionIDs){
        ids.push_back(id);
    }
    return ids;
}

std::vector<int> SOrigin::GetTracksConnectionsIDs(){
    std::vector<int> ids;
    for(const int &id: fTrackConnectionIDs){
        ids.push_back(id);
    }
    return ids;
}

void SOrigin::AddTrack(SLinearCluster track, SPoint p, double yError){
    fTrackList.push_back(track);
    fMultiplicity++;
    SPoint newPoint = SPoint( (fVertex.X()+p.X())/2., (fVertex.Y()+p.Y())/2.  );
    fVertex = newPoint;
    fYError = std::sqrt( std::pow(fYError,2) + std::pow(yError,2) );
    fNHits+=track.NHits();
}

bool SOrigin::HasTrackIndex(int ix){
    bool hasIx=false;
    for(SLinearCluster &trk: fTrackList){
        if(ix==trk.GetId()) hasIx=true;
    }
    return hasIx;
}

bool SOrigin::IsConnectedToTrackIndex(int ix){
    bool isConnected=false;
    for(const int &id:fTrackConnectionIDs){
        if(id==ix) isConnected=true;
    }
    return isConnected;
}

double SOrigin::TotalCharge(){
    double totalCharge = 0;
    for(SLinearCluster &trk: fTrackList){
        for(SHit &hit:trk.GetHits()){
            totalCharge+=hit.Integral();
        }
    }
    return totalCharge;
}

SEvent::SEvent(std::vector<SLinearCluster> tracks, std::vector<SOrigin> origins, std::vector<STriangle> angles,  std::vector<SOrigin> associatedOrigins, double hitDensity, std::vector<SHit> freeHits):
    fTrackList(tracks),
    fOriginList(origins),
    fAngleList(angles),
    fAssociatedOrigins(associatedOrigins),
    fHitDensity(hitDensity),
    fFreeHits(freeHits)
{
    // Fill the track connections
    for(SLinearCluster &track:fTrackList){
        fTrackConnections[track.GetId()] = {};
    }
    // Loop over origins
    for(SOrigin & ori:fOriginList){
        // track ids associated to the origin
        std::vector<int> trackIDs = ori.GetTracksIDs();
        for(int & id1:trackIDs){
            for(int & id2:trackIDs){
                fTrackConnections[std::abs(id1)].push_back( std::abs(id2) );
            }
        }
    }
}

int SEvent::GetNOrigins(){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(!ori.IsEdgeOrigin()) continue;
        n++;
    }
    return n;
}

int SEvent::GetNOriginsMult(int mult){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(!ori.IsEdgeOrigin()) continue;
        if(ori.Multiplicity()==mult){
            n++;
        }
    }
    return n;
}

int SEvent::GetNOriginsMultGt(int mult){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(ori.Multiplicity()>=mult){
            n++;
        }
    }
    return n;
}

int SEvent::GetNOriginsMultGt(int mult, int id1, int id2){
    int n=0;
    for(SOrigin & ori :fOriginList){
        if(ori.Multiplicity()>=mult){
            if( !(ori.HasTrackIndex(id1) || ori.HasTrackIndex(id2))){
                n++;
            }
        }
    }
    return n;
}

int SEvent::NHits(){
    int n=0;
    for(SLinearCluster & trk: fTrackList){
        n+=trk.NHits();
    }
    return n;
}

void SEvent::PrintTrackConnections(){
    for(auto & it:fTrackConnections){
        std::cout<<"Track "<<it.first<<" is connected to tracks: ";
        for(int & id:it.second){
            std::cout<<id<<" ";
        }
        std::cout<<std::endl;
    }
}

bool SEvent::IsOriginAssociatedToTrack(SOrigin ori, int trackId){
    bool isAssociated=false;
    // get the connections for trackID
    std::vector<int> connections = fTrackConnections[trackId];

    // loop over the connections
    for(int & id:ori.GetTracksConnectionsIDs()){
        for(int & id2:connections){
            if(id==id2) isAssociated=true;
        }
    }
   
    return isAssociated;
}

bool SEvent::IsTrackAssociatedToTrack(int trackId, int trackId2){
    bool isAssociated=false;
    // get the connections for trackID
    std::vector<int> connections = fTrackConnections[trackId];

    // loop over the connections
    for(int & id:connections){
        if(id==trackId2) isAssociated=true;
    }

    return isAssociated;
}


// Get origins not associated to the V+single track
SEvent SEvent::UnassociatedOrigins(STriangle triangle){

    std::vector<SOrigin> origins = GetOrigins();
    std::vector<SOrigin> notAssociatedOrigins;

    int track1ID = triangle.GetTrack1().GetId();
    int track2ID = triangle.GetTrack2().GetId();
    int mainTrackID = triangle.GetMainTrack().GetId();
    for(SOrigin &ori:origins){
        if(!IsOriginAssociatedToTrack(ori, track1ID)
            && !IsOriginAssociatedToTrack(ori, track2ID)
            && !IsOriginAssociatedToTrack(ori, mainTrackID)){
            notAssociatedOrigins.push_back(ori);
        }
    }

    SEvent notAddociatedRecoEvent({}, notAssociatedOrigins, {}, {}, 0, {});

    return notAddociatedRecoEvent;
}


void SEvent::FreeHitsAroundTriangle(STriangle triangle,
                                    int & nDirtHitsInTriangle,
                                    double & nFractionDirtHitsInTriangle,
                                    int & nDirtHitsInTriangleWires,
                                    double & nFractionDirtHitsInTriangleWires){
    // get the number of hits around the triangle
    nDirtHitsInTriangle = 0;
    nFractionDirtHitsInTriangle = 0;
    nDirtHitsInTriangleWires = 0;
    nFractionDirtHitsInTriangleWires = 0;
   
    int track1ID = triangle.GetTrack1().GetId();
    int track2ID = triangle.GetTrack2().GetId();
    int mainTrackID = triangle.GetMainTrack().GetId();
    
    // vector with the hits not assocaited t the V+main track
    std::vector<SHit> otherHits;
    std::vector<SHit> auxHits = GetFreeHits();
    otherHits.insert(otherHits.end(), auxHits.begin(), auxHits.end()); 
    for(SLinearCluster & track:GetTracks()){
        if(track.GetId()==track1ID) continue;
        if(track.GetId()==track2ID) continue;
        if(track.GetId()==mainTrackID) continue;

        // check the track is not associated to the V
        if(IsTrackAssociatedToTrack(track.GetId(), track1ID)) continue;
        if(IsTrackAssociatedToTrack(track.GetId(), track2ID)) continue;

        auxHits.clear();
        auxHits = track.GetHits();
        otherHits.insert(otherHits.end(), auxHits.begin(), auxHits.end());
    }

    
    // loop over hits, get the ones inside the triangle limits
    double triangleYLength = triangle.GetMaxY() -  triangle.GetMinY();
    for(SHit & h:otherHits){
        if(h.X()>triangle.GetMinX() && h.X()<triangle.GetMaxX() 
            && h.Y()>triangle.GetMinY() && h.Y()<triangle.GetMaxY() ){
            nDirtHitsInTriangle++;
        }

        if(h.X()>triangle.GetMinX() && h.X()<triangle.GetMaxX() 
            && h.Y()>triangle.GetMinY()-triangleYLength/2. && h.Y()<triangle.GetMaxY() + triangleYLength/2. ){
            nDirtHitsInTriangleWires++;
        }
        
    }

    nFractionDirtHitsInTriangle = 1.*nDirtHitsInTriangle/triangle.GetNHitsTriangle();
    nFractionDirtHitsInTriangleWires = 1.*nDirtHitsInTriangleWires/triangle.GetNHitsTriangle();

    return;
   
}


// Function to get free and unassociated hits in the range of the V+line
void SEvent::GetUnassociatedHits(STriangle triangle, int &nFreeHits, int &nUnassociatedHits){

    nFreeHits = 0;
    nUnassociatedHits = 0;

    // IDs of the tracks associated to the V+line
    int track1ID = triangle.GetTrack1().GetId();
    int track2ID = triangle.GetTrack2().GetId();
    int mainTrackID = triangle.GetMainTrack().GetId();

    // IDs of the tracks associated to the V+line
    std::unordered_set<int> vetoTrackIDs = {track1ID, track2ID, mainTrackID};
    std::unordered_set<int> vetoClusterIDs = {triangle.GetTrack1().GetHitCluster().GetMainClusterID(),
                                                triangle.GetTrack2().GetHitCluster().GetMainClusterID(),
                                                triangle.GetMainTrack().GetHitCluster().GetMainClusterID()};

    for(SLinearCluster & track:GetTracks()){
        if(IsTrackAssociatedToTrack(track.GetId(), track1ID)){
            vetoTrackIDs.insert(track.GetId());
            vetoClusterIDs.insert(track.GetHitCluster().GetMainClusterID());
        }
        if(IsTrackAssociatedToTrack(track.GetId(), track2ID)){
            vetoTrackIDs.insert(track.GetId());
            vetoClusterIDs.insert(track.GetHitCluster().GetMainClusterID());
        }
        if(IsTrackAssociatedToTrack(track.GetId(), mainTrackID)){
            vetoTrackIDs.insert(track.GetId());
            vetoClusterIDs.insert(track.GetHitCluster().GetMainClusterID());
        }
    }
    
    
    // Organize free hits by cluster ID
    std::unordered_set<int> freeClusterIDs;
    std::map<int, std::vector<SHit>> hitCluterIDMap;
    for(SHit & h:GetFreeHits())
        freeClusterIDs.insert(h.ClusterId());
    for(auto &id:freeClusterIDs)
        hitCluterIDMap[id] = {};
    for(SHit & h:GetFreeHits())
        hitCluterIDMap[h.ClusterId()].push_back(h);
    

    // Michel electron check
    // check the cluster doesn't start in the end of the cluster associated to the muon
    std::vector<SHit> muonTrackHits = triangle.GetMainTrack().GetHits();
    int muonClusterID = triangle.GetMainTrack().GetHitCluster().GetMainClusterID();
    if(hitCluterIDMap.find(muonClusterID)!=hitCluterIDMap.end()){
        muonTrackHits.insert(muonTrackHits.end(), hitCluterIDMap[muonClusterID].begin(), hitCluterIDMap[muonClusterID].end());
    }

    // muon end hit
    double maxD = 0;
    SHit muonEndHit=muonTrackHits.front();
    for(SHit &h:muonTrackHits){
        double d = std::abs( h.X()-triangle.GetMainVertex().X());
        if(d>maxD) {
            maxD = d;
            muonEndHit = h;
        }
    }

    // check min distance of the other cluster to the end hit (not already in the veto)
    for(auto &clusterHits:hitCluterIDMap){
        if(vetoClusterIDs.find(clusterHits.first)==vetoClusterIDs.end()){            
            
            // now check the min distance to the end hit
            double minD = 1e4;
            for(SHit &h:clusterHits.second){
                double d = std::hypot( 0.3*(h.X()-muonEndHit.X()), 0.075*(h.Y()-muonEndHit.Y()) );
                if(d<minD) minD = d;
            }

            // WARNING: HARCODED!!!!
            if(minD<2) vetoClusterIDs.insert(clusterHits.first);
        }
    }
    
    for(SLinearCluster & track:GetTracks()){
        if( vetoTrackIDs.find(track.GetId()) != vetoTrackIDs.end() ) continue;
        
        std::vector<SHit> auxHits;
        auxHits = track.GetHits();
        for(SHit & h:auxHits){
            if(vetoClusterIDs.find(h.ClusterId())==vetoClusterIDs.end()) nUnassociatedHits++;    
        }
    }

    // finally add the free hits
    for(SHit &h:GetFreeHits()){
        if( vetoClusterIDs.find( h.ClusterId() ) == vetoClusterIDs.end() )  nFreeHits++;
    }
    nUnassociatedHits+=nFreeHits;

    return;
}