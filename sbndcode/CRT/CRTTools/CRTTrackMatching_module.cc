/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTTrackMatching
/// Module Type: producer
/// File:        CRTTrackMatching_module.cc
///
/// Author:         Thomas Brooks
/// E-mail address: tbrooks@fnal.gov
///
/// Modified from CRTTrackMatching by Thomas Warburton.
/////////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/RecoUtils/RecoUtils.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

// ROOT
#include "TVector3.h"

namespace {
  // Local namespace for local functions

}

namespace sbnd {


  struct RecoCRTTrack{
    int crtID;
    int tpc;
    TVector3 start;
    TVector3 end;
    double trueTime;
    bool complete;
  };

  
  class CRTTrackMatching : public art::EDProducer {
  public:

    explicit CRTTrackMatching(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTTrackMatching(CRTTrackMatching const &) = delete;
    CRTTrackMatching(CRTTrackMatching &&) = delete;
    CRTTrackMatching & operator = (CRTTrackMatching const &) = delete; 
    CRTTrackMatching & operator = (CRTTrackMatching &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);

    // Function to transform a CRTTrack into an expected reconstructed track
    std::vector<RecoCRTTrack> CrtToRecoTrack(crt::CRTTrack, int id);

    // Function to shift CRTTrack in X and work out how much is reconstructed
    std::vector<RecoCRTTrack> CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                                 int tpc, int id, double time, bool complete);

    // Function to calculate if a CRTTrack crosses the TPC volume
    bool CrossesTPC(crt::CRTTrack track);

  private:

    // Params got from fcl file.......
    art::InputTag fTpcTrackModuleLabel; ///< name of track producer
    art::InputTag fCrtTrackModuleLabel; ///< name of crt producer
    bool          fVerbose;             ///< print info
    double        fMaxAngleDiff;        ///< max difference between CRT and TPC angles
    double        fMaxDistance;         ///< max distance between CRT and TPC start/end positions

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;              ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties; ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks; ///< pointer to detector properties provider

  }; // class CRTTrackMatching


  CRTTrackMatching::CRTTrackMatching(fhicl::ParameterSet const & p)
  // Initialize member data here, if know don't want to reconfigure on the fly
  {
    // Call appropriate produces<>() functions here.
    produces< std::vector<anab::T0> >();
    produces< art::Assns<recob::Track , anab::T0> >();
    
    // Get a pointer to the geometry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

    reconfigure(p);
  } // CRTTrackMatching()


  void CRTTrackMatching::reconfigure(fhicl::ParameterSet const & p)
  {

    fTpcTrackModuleLabel = (p.get<art::InputTag> ("TpcTrackModuleLabel"));
    fCrtTrackModuleLabel = (p.get<art::InputTag> ("CrtTrackModuleLabel")); 
    fVerbose             = (p.get<bool>          ("Verbose"));
    fMaxAngleDiff        = (p.get<double>        ("MaxAngleDiff"));
    fMaxDistance         = (p.get<double>        ("MaxDistance"));

  } // CRTTrackMatching::reconfigure()


  void CRTTrackMatching::beginJob()
  {
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;

  } // CRTTrackMatching::beginJob()


  void CRTTrackMatching::produce(art::Event & event)
  {

    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }
    
    // Create anab::T0 objects and make association with recob::Track
    std::unique_ptr< std::vector<anab::T0> > T0col( new std::vector<anab::T0>);
    std::unique_ptr< art::Assns<recob::Track, anab::T0> > Trackassn( new art::Assns<recob::Track, anab::T0>);

    // Get TPC tracks
    art::Handle< std::vector<recob::Track> > tpcTrackListHandle;
    std::vector<art::Ptr<recob::Track> > tpcTrackList;
    if (event.getByLabel(fTpcTrackModuleLabel, tpcTrackListHandle))
      art::fill_ptr_vector(tpcTrackList, tpcTrackListHandle);   

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackListHandle, event, fTpcTrackModuleLabel);

    // Get CRT tracks
    art::Handle< std::vector<crt::CRTTrack> > crtTrackListHandle;
    std::vector<art::Ptr<crt::CRTTrack> > crtTrackList;
    if (event.getByLabel(fCrtTrackModuleLabel, crtTrackListHandle))
      art::fill_ptr_vector(crtTrackList, crtTrackListHandle);

    // Validity check
    if (tpcTrackListHandle.isValid() && crtTrackListHandle.isValid() ){

      if(fVerbose) std::cout<<"Number of TPC tracks = "<<tpcTrackList.size()<<"\n"
                            <<"Number of CRT tracks = "<<crtTrackList.size()<<"\n";

      //TODO Account for crt track errors
      int crtIndex = 0;
      std::vector<RecoCRTTrack> recoCrtTracks;

      // Transform CRT tracks to expected TPC reconstructed tracks
      for (size_t crt_i = 1; crt_i < crtTrackList.size(); crt_i++){
        crt::CRTTrack crtTrack = *crtTrackList[crt_i];

        //Check that crt track crosses tpc volume, if not skip it
        if(!CrossesTPC(crtTrack)){ crtIndex++; continue; }
 
        std::vector<RecoCRTTrack> tempTracks = CrtToRecoTrack(crtTrack, crtIndex);
        recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
 
        crtIndex++;
      }
 
      // Loop over the reco crt tracks
      for (auto const& recoCrtTrack : recoCrtTracks){

        std::vector<std::pair<int, double>> crtTpcMatchCandidates;
        // Loop over the TPC tracks
        for (size_t tpc_i = 0; tpc_i < tpcTrackList.size(); tpc_i++){
          recob::Track tpcTrack = *tpcTrackList[tpc_i];

          int trackID = tpcTrack.ID();
          // If the tpcTrack has been stitched across the CPA it already has an associated t0
          std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
          int tpc = hits[0]->WireID().TPC;
 
          // Get the length, angle and start and end position of the TPC track
          TVector3 tpcStart = tpcTrack.Vertex();
          TVector3 tpcEnd = tpcTrack.End();
          double tpcTheta = (tpcStart - tpcEnd).Theta();
          double tpcPhi = (tpcStart - tpcEnd).Phi();
 
          // Get the length, angle and start and end position of the TPC track
          TVector3 crtStart = recoCrtTrack.start;
          TVector3 crtEnd = recoCrtTrack.end;
          double crtTheta = (crtStart - crtEnd).Theta();
          double crtPhi = (crtStart - crtEnd).Phi();
 
          // Find the difference with the CRT track
          double dDist1 = (crtStart-tpcStart).Mag();
          double dDist2 = (crtEnd-tpcEnd).Mag();
          if(std::max((crtStart-tpcStart).Mag(), (crtEnd-tpcEnd).Mag()) > 
             std::max((crtStart-tpcEnd).Mag(), (crtEnd-tpcStart).Mag())){
            crtTheta = (crtEnd - crtStart).Theta();
            crtPhi = (crtEnd - crtStart).Phi();
            dDist1 = (crtEnd-tpcStart).Mag();
            dDist2 = (crtStart-tpcEnd).Mag();
          }
          double dTheta = atan2(sin(tpcTheta - crtTheta), cos(tpcTheta - crtTheta));
          double dPhi = atan2(sin(tpcPhi - crtPhi), cos(tpcPhi - crtPhi));
 
          // Do the actual matching
          if(std::abs(dTheta) < fMaxAngleDiff && std::abs(dPhi) < fMaxAngleDiff && 
             tpc == recoCrtTrack.tpc && (dDist1<fMaxDistance||dDist2<fMaxDistance)){
            crtTpcMatchCandidates.push_back(std::make_pair(trackID, std::abs(dTheta)));
          }
             
        }
        // Choose the track which matches the closest
        int matchedTrackID = -99999;
        if(crtTpcMatchCandidates.size() > 0){
          std::sort(crtTpcMatchCandidates.begin(), crtTpcMatchCandidates.end(), [](auto& left, auto& right){
                    return left.second < right.second;});
          matchedTrackID = crtTpcMatchCandidates[0].first;
          if(fVerbose) std::cout<<"Matched time "<<recoCrtTrack.trueTime<<" ticks to track "<<matchedTrackID<<"\n";
        }
        if(matchedTrackID != -99999){
          T0col->push_back(anab::T0(fDetectorClocks->TPCTick2Time(recoCrtTrack.trueTime)*1e3, 0, 
                                    matchedTrackID, (*T0col).size(), crtTpcMatchCandidates[0].second)); 
          util::CreateAssn(*this, event, *T0col, tpcTrackList[matchedTrackID], *Trackassn);
        }
 
      }

    } // Validity check
   
    event.put(std::move(T0col));
    event.put(std::move(Trackassn));
    
  } // CRTTrackMatching::produce()


  void CRTTrackMatching::endJob()
  {

  } // CRTTrackMatching::endJob()

  
  // Function to transform a CRTTrack into an expected reconstructed track
  std::vector<RecoCRTTrack> CRTTrackMatching::CrtToRecoTrack(crt::CRTTrack track, int id){

    // TODO: could be more detector agnositic (pos of cpa, drift directions, max length of track)
    std::vector<RecoCRTTrack> recoCrtTracks;
    // Get the time of the track
    double crtTimeTicks = fDetectorClocks->TPCG4Time2Tick((double)(int)track.ts1_ns); // ns -> tick
    // Convert time into a x shift
    double xShift = fDetectorClocks->TPCTick2Time(crtTimeTicks) * fDetectorProperties->DriftVelocity();

    // Shift track, remembering to take into account the tpc, if the track crosses the cpa and 
    //the size of the readout window
    TVector3 crtStart(track.x1_pos, track.y1_pos, track.z1_pos);
    TVector3 crtEnd(track.x2_pos, track.y2_pos, track.z2_pos);
    if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
    if(!track.complete){
      // Find point where track crosses bottom plane
      TVector3 diff = (crtEnd - crtStart).Unit();
      TVector3 newEnd = crtStart + 5000*diff; 
      crtEnd = newEnd;
    }

    TVector3 cpaCrossStart, cpaCrossEnd;
    // Calculate the expected reconstructed length, angle and start and end position
    std::vector<int> crtTpc;
    if(crtStart.X() < 0. && crtEnd.X() < 0.){

      // Track in TPC 0
      std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 0, 
                                                                id, crtTimeTicks, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
    }
    else if(crtStart.X() > 0. && crtEnd.X() > 0.){

      // Track in TPC 1
      std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 1, 
                                                                id, crtTimeTicks, track.complete);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
    }
    else {
      // Track in both TPCs and will be split
      TVector3 direction = crtStart - crtEnd;
      double step = (0. - crtStart.X())/direction.X(); 
      TVector3 cpaCross(0., crtStart.Y() + step*direction.Y(), crtStart.Z() + step*direction.Z());
      cpaCrossStart = cpaCross;
      cpaCrossEnd = cpaCross;

      if(crtStart.X() < 0.){ 
        std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 0, 
                                                                   id, crtTimeTicks, track.complete);
        recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

        std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 1, 
                                                                   id, crtTimeTicks, track.complete);
        recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
      }
      else {
        std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 0, 
                                                                   id, crtTimeTicks, track.complete);
        recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

        std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 1, 
                                                                   id, crtTimeTicks, track.complete);
        recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
      }

    }

    return recoCrtTracks;

  } // CRTTrackMatching::CRTToRecoTrack()


  // Function to shift CRTTrack in X and work out how much is reconstructed
  std::vector<RecoCRTTrack> CRTTrackMatching::CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                                                 int tpc, int id, double time, bool complete){

    std::vector<RecoCRTTrack> recoCrtTracks;

    // Get the true entry and exit points in the TPC
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();

    // Get track info
    TVector3 diff = end - start;
    TVector3 startTPC, endTPC;
    bool first = true;

    // Loop over trajectory points
    int npts = 1000;
    for (int traj_i = 0; traj_i <= npts; traj_i++){
      TVector3 trajPoint = start + ((traj_i)/((double)npts))*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
          && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
        first = false;
        startTPC = trajPoint;
      }
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
          && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
        endTPC = trajPoint;
      }
    }

    // Don't shift if not inside TPC
    if(startTPC == endTPC){
      return recoCrtTracks;
    }

    // Shift in x depending on TPC
    if(tpc == 1){
      // Track in TPC 1
      startTPC[0] -= shift;
      endTPC[0] -= shift;
    }
    else if(tpc == 0){
      // Track in TPC 0
      startTPC[0] += shift;
      endTPC[0] += shift;
    }
    
    double readoutWindow  = (double)fDetectorProperties->ReadOutWindowSize();
    double driftTimeTicks = fDetectorClocks->Time2Tick((2.*fGeometryService->DetHalfWidth())/fDetectorProperties->DriftVelocity());
    double deltaX = fDetectorClocks->TPCTick2Time(readoutWindow - driftTimeTicks) * fDetectorProperties->DriftVelocity();

    if(tpc == 0) xmax = deltaX;
    if(tpc == 1) xmin = -deltaX;

    // Get track info
    TVector3 diffTPC = endTPC - startTPC;
    TVector3 startCut, endCut;
    first = true;
    // Loop over trajectory points
    for (int traj_i = 0; traj_i <= npts; traj_i++){
      TVector3 trajPoint = startTPC + ((traj_i)/((double)npts))*diffTPC;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
          && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
        first = false;
        startCut = trajPoint;
      }
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
          && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
        endCut = trajPoint;
      }
    }

    if(!(startCut.X() == endCut.X())){
      RecoCRTTrack recoCrtTrack = {id, tpc, startCut, endCut, time, complete};
      recoCrtTracks.push_back(recoCrtTrack);
    }

    return recoCrtTracks;

  } // CRTTrackMatching::CreateRecoCRTTrack()


  // Function to calculate if a CRTTrack crosses the TPC volume
  bool CRTTrackMatching::CrossesTPC(crt::CRTTrack track){

    // Check if particle enters the TPC
    bool enters = false;
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();
    if(track.complete){
      // Get track info
      TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
      TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
      TVector3 diff = end - start;
      // Loop over trajectory points
      int npts = 100;
      for (int traj_i = 0; traj_i < npts; traj_i++){
        TVector3 trajPoint = start + ((traj_i+1)/100.)*diff;
        // Check if point is within reconstructable volume
        if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
          enters = true;
        }
      }
    }
    // If track just between top two planes
    else{
      //
      TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
      TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
      if(start.Y() < end.Y()) std::swap(start, end);
      TVector3 diff = (end - start).Unit();
      int npts = 100;
      for (int traj_i = 0; traj_i < npts; traj_i++){
        // TODO: length of track needs to be more detector agnostic
        TVector3 trajPoint = start + (100.*(traj_i+1))*diff;
        // Check if point is within reconstructable volume
        if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
          enters = true;
        }
      }
    }

    return enters;

  } // CRTTrackMatching::CrossesTPC()


  DEFINE_ART_MODULE(CRTTrackMatching)

} // sbnd namespace

namespace {

}
